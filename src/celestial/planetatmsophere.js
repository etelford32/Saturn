/**
 * PlanetAtmosphere.js - Saturn's atmospheric dynamics and weather systems
 * Implements shallow water equations, jet streams, and the famous hexagonal polar vortex
 */

import * as THREE from 'three';

export class PlanetAtmosphere {
    constructor(planet, config = {}) {
        this.planet = planet;
        
        // Grid resolution for atmospheric simulation
        this.gridSize = config.gridSize || 128;
        this.layers = config.layers || 3;
        
        // Physical constants
        this.constants = {
            gasConstant: 3777, // J/kg⋅K for H2-He atmosphere
            adiabaticIndex: 1.4,
            coriolisParameter: 2 * (2 * Math.PI / planet.rotationPeriod), // 2Ω
            planetaryRadius: planet.radius * 1000, // meters
            gravitationalAccel: planet.getSurfaceGravity()
        };
        
        // Atmospheric grid (spherical coordinates)
        this.grid = {
            longitude: [], // φ (0 to 2π)
            latitude: [],  // θ (0 to π)
            pressure: [],  // p levels
            deltaLon: 2 * Math.PI / this.gridSize,
            deltaLat: Math.PI / this.gridSize,
            deltaPressure: []
        };
        
        // State variables on grid
        this.state = {
            // Velocity components (m/s)
            u: [], // Zonal (eastward) velocity
            v: [], // Meridional (northward) velocity
            w: [], // Vertical velocity
            
            // Thermodynamic variables
            temperature: [], // K
            pressure: [],    // Pa
            density: [],     // kg/m³
            
            // Heights
            geopotentialHeight: [], // m²/s²
            
            // Vorticity and divergence
            vorticity: [],
            divergence: [],
            
            // Potential vorticity (conserved quantity)
            potentialVorticity: []
        };
        
        // Jet stream configuration
        this.jetStreams = [
            { latitude: 80, speed: 50, width: 10 },    // North polar
            { latitude: 60, speed: -100, width: 8 },   // North temperate
            { latitude: 40, speed: 150, width: 12 },   // North subtropical
            { latitude: 20, speed: -80, width: 8 },    // North tropical
            { latitude: 0, speed: 500, width: 15 },    // Equatorial superrotation
            { latitude: -20, speed: -80, width: 8 },   // South tropical
            { latitude: -40, speed: 150, width: 12 },  // South subtropical
            { latitude: -60, speed: -100, width: 8 },  // South temperate
            { latitude: -80, speed: 50, width: 10 }    // South polar
        ];
        
        // Hexagonal wave at north pole
        this.hexagon = {
            latitude: 78, // degrees north
            waveNumber: 6, // hexagonal pattern
            amplitude: 50, // m/s velocity amplitude
            phase: 0,      // current phase
            frequency: 1 / (10.7 * 3600), // One rotation period
            meanFlow: 100  // Background westerly flow (m/s)
        };
        
        // Forcing terms
        this.forcing = {
            // Solar heating
            solarHeating: {
                solarConstant: 14.9, // W/m² at Saturn
                albedo: planet.albedo,
                absorptionEfficiency: 0.95
            },
            
            // Internal heating
            internalHeating: {
                flux: 2.01 * 14.9, // W/m² (Saturn radiates 2.01x more than it receives)
                distribution: 'uniform' // Simplified
            },
            
            // Radiative cooling
            radiativeCooling: {
                timeConstant: 86400 * 100, // 100 days (approximation)
                emissivity: 0.95
            },
            
            // Frictional damping
            friction: {
                coefficient: 1e-6, // s⁻¹
                heightScale: 10000 // m (scale height for friction)
            }
        };
        
        // Storm systems
        this.storms = [];
        this.stormGeneration = {
            probability: 1e-6, // per grid point per time step
            lifetimeRange: [3600 * 24, 3600 * 24 * 30], // 1-30 days
            intensityRange: [10, 100] // m/s maximum winds
        };
        
        // Seasonal effects (Saturn's year is ~29.5 Earth years)
        this.season = {
            currentTime: 0, // seconds since vernal equinox
            orbitalPeriod: 29.5 * 365.25 * 24 * 3600, // seconds
            axialTilt: planet.axialTilt,
            currentSolarDeclination: 0
        };
        
        // Initialize atmospheric model
        this.initializeAtmosphere();
    }
    
    /**
     * Initialize the atmospheric grid and state
     */
    initializeAtmosphere() {
        this.setupGrid();
        this.initializeState();
        this.setupJetStreams();
        this.initializeHexagon();
        this.setupTemperatureProfile();
    }
    
    /**
     * Setup the spherical grid
     */
    setupGrid() {
        // Longitude grid (0 to 2π)
        for (let i = 0; i < this.gridSize; i++) {
            this.grid.longitude.push(i * this.grid.deltaLon);
        }
        
        // Latitude grid (0 to π, but we'll use -π/2 to π/2 for convenience)
        for (let j = 0; j < this.gridSize; j++) {
            this.grid.latitude.push(-Math.PI/2 + j * this.grid.deltaLat);
        }
        
        // Pressure levels (log spacing)
        const pTop = 1000; // Pa (top of model)
        const pBottom = 1e6; // Pa (bottom of model)
        for (let k = 0; k < this.layers; k++) {
            const pressure = pTop * Math.pow(pBottom/pTop, k/(this.layers-1));
            this.grid.pressure.push(pressure);
        }
    }
    
    /**
     * Initialize atmospheric state variables
     */
    initializeState() {
        const numPoints = this.gridSize * this.gridSize * this.layers;
        
        // Initialize arrays
        this.state.u = new Float32Array(numPoints);
        this.state.v = new Float32Array(numPoints);
        this.state.w = new Float32Array(numPoints);
        this.state.temperature = new Float32Array(numPoints);
        this.state.pressure = new Float32Array(numPoints);
        this.state.density = new Float32Array(numPoints);
        this.state.geopotentialHeight = new Float32Array(numPoints);
        this.state.vorticity = new Float32Array(numPoints);
        this.state.divergence = new Float32Array(numPoints);
        this.state.potentialVorticity = new Float32Array(numPoints);
        
        // Set initial conditions
        for (let k = 0; k < this.layers; k++) {
            for (let j = 0; j < this.gridSize; j++) {
                for (let i = 0; i < this.gridSize; i++) {
                    const index = this.getGridIndex(i, j, k);
                    const lat = this.grid.latitude[j];
                    const pressure = this.grid.pressure[k];
                    
                    // Initial temperature profile
                    this.state.temperature[index] = this.getInitialTemperature(lat, pressure);
                    
                    // Hydrostatic pressure
                    this.state.pressure[index] = pressure;
                    
                    // Ideal gas law for density
                    this.state.density[index] = pressure / (this.constants.gasConstant * this.state.temperature[index]);
                    
                    // Initial winds (small perturbations)
                    this.state.u[index] = (Math.random() - 0.5) * 2; // ±1 m/s
                    this.state.v[index] = (Math.random() - 0.5) * 2;
                    this.state.w[index] = 0;
                }
            }
        }
    }
    
    /**
     * Get initial temperature based on latitude and pressure
     */
    getInitialTemperature(latitude, pressure) {
        const cosLat = Math.cos(latitude);
        
        // Basic temperature profile
        let temperature = 95; // Base temperature (K)
        
        // Latitudinal variation
        temperature += 20 * cosLat * cosLat; // Warmer at equator
        
        // Pressure variation (lapse rate)
        const refPressure = 1e5; // Pa
        const lapseRate = 2e-3; // K/Pa
        temperature += lapseRate * (pressure - refPressure);
        
        return Math.max(temperature, 50); // Minimum 50K
    }
    
    /**
     * Setup initial jet stream pattern
     */
    setupJetStreams() {
        this.jetStreams.forEach(jet => {
            const jetLat = THREE.MathUtils.degToRad(jet.latitude);
            const jetWidth = THREE.MathUtils.degToRad(jet.width);
            
            for (let k = 0; k < this.layers; k++) {
                for (let j = 0; j < this.gridSize; j++) {
                    for (let i = 0; i < this.gridSize; i++) {
                        const index = this.getGridIndex(i, j, k);
                        const lat = this.grid.latitude[j];
                        
                        // Gaussian profile for jet
                        const latDiff = lat - jetLat;
                        const jetProfile = Math.exp(-0.5 * Math.pow(latDiff / jetWidth, 2));
                        
                        // Add jet contribution to zonal wind
                        this.state.u[index] += jet.speed * jetProfile;
                    }
                }
            }
        });
    }
    
    /**
     * Initialize the hexagonal polar vortex
     */
    initializeHexagon() {
        const hexLat = THREE.MathUtils.degToRad(this.hexagon.latitude);
        const hexWidth = THREE.MathUtils.degToRad(10); // 10 degree width
        
        for (let k = 0; k < this.layers; k++) {
            for (let j = 0; j < this.gridSize; j++) {
                for (let i = 0; i < this.gridSize; i++) {
                    const index = this.getGridIndex(i, j, k);
                    const lat = this.grid.latitude[j];
                    const lon = this.grid.longitude[i];
                    
                    // Only apply near hexagon latitude
                    const latProfile = Math.exp(-0.5 * Math.pow((lat - hexLat) / hexWidth, 2));
                    
                    if (latProfile > 0.1) {
                        // Hexagonal wave perturbation
                        const hexWave = this.hexagon.amplitude * Math.cos(this.hexagon.waveNumber * lon + this.hexagon.phase);
                        
                        // Add to zonal wind
                        this.state.u[index] += (this.hexagon.meanFlow + hexWave) * latProfile;
                    }
                }
            }
        }
    }
    
    /**
     * Setup realistic temperature profile
     */
    setupTemperatureProfile() {
        for (let k = 0; k < this.layers; k++) {
            for (let j = 0; j < this.gridSize; j++) {
                for (let i = 0; i < this.gridSize; i++) {
                    const index = this.getGridIndex(i, j, k);
                    const lat = this.grid.latitude[j];
                    const pressure = this.grid.pressure[k];
                    
                    // Update temperature with more realistic profile
                    this.state.temperature[index] = this.calculateRealisticTemperature(lat, pressure);
                    
                    // Update density
                    this.state.density[index] = this.state.pressure[index] / 
                        (this.constants.gasConstant * this.state.temperature[index]);
                }
            }
        }
    }
    
    /**
     * Calculate realistic temperature profile
     */
    calculateRealisticTemperature(latitude, pressure) {
        const absLat = Math.abs(latitude);
        
        // Base temperature from observations
        let temp = 95; // Tropospheric temperature
        
        // Latitudinal gradient
        temp += 30 * Math.cos(latitude) * Math.cos(latitude); // Equator warmer
        
        // Polar cooling
        if (absLat > Math.PI/3) { // Beyond 60 degrees
            temp -= 20 * Math.pow((absLat - Math.PI/3) / (Math.PI/2 - Math.PI/3), 2);
        }
        
        // Pressure dependence (stratosphere heating)
        if (pressure < 1e4) { // Above troposphere
            const stratosphericHeating = 65 * Math.exp(-(pressure - 1e3) / 5e3);
            temp += stratosphericHeating;
        }
        
        return Math.max(temp, 40); // Minimum temperature
    }
    
    /**
     * Get grid index from coordinates
     */
    getGridIndex(i, j, k) {
        return k * this.gridSize * this.gridSize + j * this.gridSize + i;
    }
    
    /**
     * Get grid coordinates from index
     */
    getGridCoords(index) {
        const k = Math.floor(index / (this.gridSize * this.gridSize));
        const remainder = index % (this.gridSize * this.gridSize);
        const j = Math.floor(remainder / this.gridSize);
        const i = remainder % this.gridSize;
        return { i, j, k };
    }
    
    /**
     * Calculate Coriolis parameter at latitude
     */
    getCoriolisParameter(latitude) {
        return this.constants.coriolisParameter * Math.sin(latitude);
    }
    
    /**
     * Calculate pressure gradient force
     */
    calculatePressureGradientForce(i, j, k) {
        const index = this.getGridIndex(i, j, k);
        const density = this.state.density[index];
        
        // Pressure gradients
        const dpdx = this.calculateLongitudinalGradient(this.state.pressure, i, j, k);
        const dpdy = this.calculateLatitudinalGradient(this.state.pressure, i, j, k);
        
        // Convert to acceleration (m/s²)
        const lat = this.grid.latitude[j];
        const radius = this.constants.planetaryRadius;
        
        const ax = -dpdx / (density * radius * Math.cos(lat));
        const ay = -dpdy / (density * radius);
        
        return { ax, ay };
    }
    
    /**
     * Calculate longitudinal gradient
     */
    calculateLongitudinalGradient(field, i, j, k) {
        const im1 = (i - 1 + this.gridSize) % this.gridSize;
        const ip1 = (i + 1) % this.gridSize;
        
        const valueLeft = field[this.getGridIndex(im1, j, k)];
        const valueRight = field[this.getGridIndex(ip1, j, k)];
        
        return (valueRight - valueLeft) / (2 * this.grid.deltaLon);
    }
    
    /**
     * Calculate latitudinal gradient
     */
    calculateLatitudinalGradient(field, i, j, k) {
        if (j === 0 || j === this.gridSize - 1) return 0; // Poles
        
        const valueNorth = field[this.getGridIndex(i, j + 1, k)];
        const valueSouth = field[this.getGridIndex(i, j - 1, k)];
        
        return (valueNorth - valueSouth) / (2 * this.grid.deltaLat);
    }
    
    /**
     * Calculate vorticity
     */
    calculateVorticity(i, j, k) {
        // ζ = ∂v/∂x - ∂u/∂y
        const dvdx = this.calculateLongitudinalGradient(this.state.v, i, j, k);
        const dudy = this.calculateLatitudinalGradient(this.state.u, i, j, k);
        
        const lat = this.grid.latitude[j];
        const radius = this.constants.planetaryRadius;
        
        return (dvdx / (radius * Math.cos(lat)) - dudy / radius);
    }
    
    /**
     * Calculate divergence
     */
    calculateDivergence(i, j, k) {
        // ∇·v = ∂u/∂x + ∂v/∂y
        const dudx = this.calculateLongitudinalGradient(this.state.u, i, j, k);
        const dvdy = this.calculateLatitudinalGradient(this.state.v, i, j, k);
        
        const lat = this.grid.latitude[j];
        const radius = this.constants.planetaryRadius;
        
        return (dudx / (radius * Math.cos(lat)) + dvdy / radius);
    }
    
    /**
     * Apply solar heating
     */
    applySolarHeating(i, j, k, deltaTime) {
        const index = this.getGridIndex(i, j, k);
        const lat = this.grid.latitude[j];
        
        // Solar zenith angle (simplified)
        const solarDeclination = this.season.currentSolarDeclination;
        const cosZenith = Math.sin(lat) * Math.sin(solarDeclination) + 
                         Math.cos(lat) * Math.cos(solarDeclination);
        
        if (cosZenith > 0) {
            const solarFlux = this.forcing.solarHeating.solarConstant * cosZenith;
            const absorbed = solarFlux * (1 - this.forcing.solarHeating.albedo) * 
                           this.forcing.solarHeating.absorptionEfficiency;
            
            // Convert to temperature change
            const heatCapacity = this.constants.gasConstant / (this.constants.adiabaticIndex - 1);
            const deltaT = absorbed * deltaTime / (this.state.density[index] * heatCapacity);
            
            this.state.temperature[index] += deltaT;
        }
    }
    
    /**
     * Apply internal heating
     */
    applyInternalHeating(i, j, k, deltaTime) {
        const index = this.getGridIndex(i, j, k);
        
        const internalFlux = this.forcing.internalHeating.flux;
        const heatCapacity = this.constants.gasConstant / (this.constants.adiabaticIndex - 1);
        const deltaT = internalFlux * deltaTime / (this.state.density[index] * heatCapacity);
        
        this.state.temperature[index] += deltaT;
    }
    
    /**
     * Apply radiative cooling
     */
    applyRadiativeCooling(i, j, k, deltaTime) {
        const index = this.getGridIndex(i, j, k);
        const currentTemp = this.state.temperature[index];
        
        // Simple Newtonian cooling toward equilibrium
        const equilibriumTemp = this.getEquilibriumTemperature(this.grid.latitude[j]);
        const coolingRate = (equilibriumTemp - currentTemp) / this.forcing.radiativeCooling.timeConstant;
        
        this.state.temperature[index] += coolingRate * deltaTime;
    }
    
    /**
     * Get equilibrium temperature for latitude
     */
    getEquilibriumTemperature(latitude) {
        const basTemp = 95; // K
        const latitudeEffect = 30 * Math.cos(latitude) * Math.cos(latitude);
        return basTemp + latitudeEffect;
    }
    
    /**
     * Apply friction damping
     */
    applyFriction(i, j, k, deltaTime) {
        const index = this.getGridIndex(i, j, k);
        const friction = this.forcing.friction.coefficient;
        
        // Exponential damping
        const dampingFactor = Math.exp(-friction * deltaTime);
        this.state.u[index] *= dampingFactor;
        this.state.v[index] *= dampingFactor;
    }
    
    /**
     * Update hexagonal polar vortex
     */
    updateHexagon(deltaTime) {
        // Advance phase
        this.hexagon.phase += this.hexagon.frequency * 2 * Math.PI * deltaTime;
        
        // Update hexagonal pattern in grid
        const hexLat = THREE.MathUtils.degToRad(this.hexagon.latitude);
        const hexWidth = THREE.MathUtils.degToRad(10);
        
        for (let k = 0; k < this.layers; k++) {
            for (let j = 0; j < this.gridSize; j++) {
                const lat = this.grid.latitude[j];
                const latProfile = Math.exp(-0.5 * Math.pow((lat - hexLat) / hexWidth, 2));
                
                if (latProfile > 0.1) {
                    for (let i = 0; i < this.gridSize; i++) {
                        const index = this.getGridIndex(i, j, k);
                        const lon = this.grid.longitude[i];
                        
                        // Remove old hexagonal contribution (approximate)
                        const oldHex = this.hexagon.amplitude * Math.cos(this.hexagon.waveNumber * lon + this.hexagon.phase - this.hexagon.frequency * 2 * Math.PI * deltaTime);
                        this.state.u[index] -= oldHex * latProfile;
                        
                        // Add new hexagonal contribution
                        const newHex = this.hexagon.amplitude * Math.cos(this.hexagon.waveNumber * lon + this.hexagon.phase);
                        this.state.u[index] += newHex * latProfile;
                    }
                }
            }
        }
    }
    
    /**
     * Generate random storms
     */
    generateStorms(deltaTime) {
        const stormProb = this.stormGeneration.probability * deltaTime;
        
        for (let k = 0; k < this.layers; k++) {
            for (let j = 0; j < this.gridSize; j++) {
                for (let i = 0; i < this.gridSize; i++) {
                    if (Math.random() < stormProb) {
                        const lat = this.grid.latitude[j];
                        const lon = this.grid.longitude[i];
                        
                        // Don't generate storms at poles
                        if (Math.abs(lat) < Math.PI/2 - 0.1) {
                            this.createStorm(lon, lat, k);
                        }
                    }
                }
            }
        }
    }
    
    /**
     * Create a storm system
     */
    createStorm(longitude, latitude, level) {
        const storm = {
            id: Date.now() + Math.random(),
            longitude: longitude,
            latitude: latitude,
            level: level,
            intensity: this.stormGeneration.intensityRange[0] + 
                      Math.random() * (this.stormGeneration.intensityRange[1] - this.stormGeneration.intensityRange[0]),
            lifetime: this.stormGeneration.lifetimeRange[0] + 
                     Math.random() * (this.stormGeneration.lifetimeRange[1] - this.stormGeneration.lifetimeRange[0]),
            age: 0,
            size: 5 + Math.random() * 15 // Grid points radius
        };
        
        this.storms.push(storm);
        this.applyStormToGrid(storm);
    }
    
    /**
     * Apply storm effects to atmospheric grid
     */
    applyStormToGrid(storm) {
        const centerI = Math.round(storm.longitude / this.grid.deltaLon);
        const centerJ = Math.round((storm.latitude + Math.PI/2) / this.grid.deltaLat);
        const size = storm.size;
        
        for (let di = -size; di <= size; di++) {
            for (let dj = -size; dj <= size; dj++) {
                const i = (centerI + di + this.gridSize) % this.gridSize;
                const j = Math.max(0, Math.min(this.gridSize - 1, centerJ + dj));
                
                const distance = Math.sqrt(di*di + dj*dj);
                if (distance <= size) {
                    const index = this.getGridIndex(i, j, storm.level);
                    const influence = Math.exp(-distance / (size/2));
                    
                    // Add cyclonic circulation
                    const angle = Math.atan2(dj, di);
                    const stormU = -storm.intensity * Math.sin(angle) * influence;
                    const stormV = storm.intensity * Math.cos(angle) * influence;
                    
                    this.state.u[index] += stormU;
                    this.state.v[index] += stormV;
                    
                    // Pressure depression
                    this.state.pressure[index] -= 1000 * influence; // Pa
                }
            }
        }
    }
    
    /**
     * Update storm systems
     */
    updateStorms(deltaTime) {
        this.storms = this.storms.filter(storm => {
            storm.age += deltaTime;
            
            // Move storm (simple eastward drift)
            storm.longitude += deltaTime * 50 / this.constants.planetaryRadius; // ~50 m/s drift
            storm.longitude %= (2 * Math.PI);
            
            // Decay intensity over time
            const decayRate = 1 / storm.lifetime;
            storm.intensity *= Math.exp(-decayRate * deltaTime);
            
            // Remove if too old or weak
            return storm.age < storm.lifetime && storm.intensity > 1;
        });
    }
    
    /**
     * Solve shallow water equations (simplified)
     */
    solveShallowWater(deltaTime) {
        // This is a simplified version - real implementation would use proper numerical schemes
        const newU = new Float32Array(this.state.u);
        const newV = new Float32Array(this.state.v);
        
        for (let k = 0; k < this.layers; k++) {
            for (let j = 1; j < this.gridSize - 1; j++) { // Skip poles
                for (let i = 0; i < this.gridSize; i++) {
                    const index = this.getGridIndex(i, j, k);
                    const lat = this.grid.latitude[j];
                    
                    // Current state
                    const u = this.state.u[index];
                    const v = this.state.v[index];
                    
                    // Coriolis parameter
                    const f = this.getCoriolisParameter(lat);
                    
                    // Pressure gradient force
                    const pgf = this.calculatePressureGradientForce(i, j, k);
                    
                    // Coriolis force
                    const coriolisU = f * v;
                    const coriolisV = -f * u;
                    
                    // Update velocities
                    newU[index] = u + deltaTime * (pgf.ax + coriolisU);
                    newV[index] = v + deltaTime * (pgf.ay + coriolisV);
                }
            }
        }
        
        // Update state
        this.state.u = newU;
        this.state.v = newV;
    }
    
    /**
     * Update seasonal effects
     */
    updateSeason(deltaTime) {
        this.season.currentTime += deltaTime;
        this.season.currentTime %= this.season.orbitalPeriod;
        
        // Calculate solar declination
        const orbitalPhase = 2 * Math.PI * this.season.currentTime / this.season.orbitalPeriod;
        this.season.currentSolarDeclination = this.season.axialTilt * Math.sin(orbitalPhase);
    }
    
    /**
     * Main atmospheric update
     */
    update(deltaTime) {
        // Update seasonal effects
        this.updateSeason(deltaTime);
        
        // Apply heating and cooling
        for (let k = 0; k < this.layers; k++) {
            for (let j = 0; j < this.gridSize; j++) {
                for (let i = 0; i < this.gridSize; i++) {
                    this.applySolarHeating(i, j, k, deltaTime);
                    this.applyInternalHeating(i, j, k, deltaTime);
                    this.applyRadiativeCooling(i, j, k, deltaTime);
                    this.applyFriction(i, j, k, deltaTime);
                }
            }
        }
        
        // Update dynamical systems
        this.solveShallowWater(deltaTime);
        this.updateHexagon(deltaTime);
        
        // Update storms
        this.generateStorms(deltaTime);
        this.updateStorms(deltaTime);
        
        // Calculate derived quantities
        this.calculateVorticityAndDivergence();
    }
    
    /**
     * Calculate vorticity and divergence fields
     */
    calculateVorticityAndDivergence() {
        for (let k = 0; k < this.layers; k++) {
            for (let j = 1; j < this.gridSize - 1; j++) {
                for (let i = 0; i < this.gridSize; i++) {
                    const index = this.getGridIndex(i, j, k);
                    this.state.vorticity[index] = this.calculateVorticity(i, j, k);
                    this.state.divergence[index] = this.calculateDivergence(i, j, k);
                }
            }
        }
    }
    
    /**
     * Get atmospheric state at position
     */
    getStateAtPosition(longitude, latitude, pressure) {
        // Find grid indices
        const i = Math.floor(longitude / this.grid.deltaLon);
        const j = Math.floor((latitude + Math.PI/2) / this.grid.deltaLat);
        const k = this.findPressureLevel(pressure);
        
        if (i < 0 || i >= this.gridSize || j < 0 || j >= this.gridSize || k < 0 || k >= this.layers) {
            return null;
        }
        
        const index = this.getGridIndex(i, j, k);
        
        return {
            u: this.state.u[index],
            v: this.state.v[index],
            w: this.state.w[index],
            temperature: this.state.temperature[index],
            pressure: this.state.pressure[index],
            density: this.state.density[index],
            vorticity: this.state.vorticity[index],
            divergence: this.state.divergence[index]
        };
    }
    
    /**
     * Find pressure level index
     */
    findPressureLevel(pressure) {
        for (let k = 0; k < this.layers - 1; k++) {
            if (pressure >= this.grid.pressure[k] && pressure <= this.grid.pressure[k + 1]) {
                return k;
            }
        }
        return this.layers - 1;
    }
    
    /**
     * Get atmospheric information
     */
    getAtmosphereInfo() {
        return {
            gridSize: this.gridSize,
            layers: this.layers,
            jetStreams: this.jetStreams,
            hexagon: this.hexagon,
            storms: this.storms.length,
            season: this.season
        };
    }
    
    /**
     * Export atmospheric state
     */
    exportState() {
        return {
            state: {
                u: Array.from(this.state.u),
                v: Array.from(this.state.v),
                temperature: Array.from(this.state.temperature),
                pressure: Array.from(this.state.pressure)
            },
            storms: this.storms,
            season: this.season,
            hexagon: this.hexagon
        };
    }
    
    /**
     * Import atmospheric state
     */
    importState(data) {
        if (data.state) {
            this.state.u = new Float32Array(data.state.u);
            this.state.v = new Float32Array(data.state.v);
            this.state.temperature = new Float32Array(data.state.temperature);
            this.state.pressure = new Float32Array(data.state.pressure);
        }
        
        if (data.storms) this.storms = data.storms;
        if (data.season) this.season = { ...this.season, ...data.season };
        if (data.hexagon) this.hexagon = { ...this.hexagon, ...data.hexagon };
    }
}
