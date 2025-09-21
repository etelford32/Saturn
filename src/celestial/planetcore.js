/**
 * PlanetCore.js - Saturn's internal structure and thermal dynamics
 * Models the complex interior structure, heat generation, and energy transport
 */

import * as THREE from 'three';

export class PlanetCore {
    constructor(planet, config = {}) {
        this.planet = planet;
        
        // Core structure layers (from center outward)
        this.layers = {
            // Inner rocky/ice core
            innerCore: {
                radius: 0.2 * planet.radius, // ~12,000 km
                composition: {
                    rock: 0.7,    // Silicates, metals
                    ice: 0.25,    // Water ice
                    organics: 0.05 // Carbon compounds
                },
                temperature: 12000, // K (very hot!)
                pressure: 8e11,     // Pa (800 GPa)
                density: 19000,     // kg/m³
                state: 'solid'
            },
            
            // Metallic hydrogen layer
            metallicHydrogen: {
                radius: 0.6 * planet.radius, // ~35,000 km
                composition: {
                    metallicH2: 0.95,
                    metallicHe: 0.05
                },
                temperature: 8000,  // K
                pressure: 4e11,     // Pa (400 GPa)
                density: 1100,      // kg/m³
                state: 'metallic'
            },
            
            // Molecular hydrogen layer
            molecularHydrogen: {
                radius: 0.85 * planet.radius, // ~49,000 km
                composition: {
                    H2: 0.89,
                    He: 0.105,
                    H2O: 0.003,
                    CH4: 0.002
                },
                temperature: 2000,  // K
                pressure: 1e9,      // Pa (1 GPa)
                density: 300,       // kg/m³
                state: 'liquid'
            },
            
            // Outer atmosphere
            atmosphere: {
                radius: planet.radius,
                composition: planet.atmosphere.composition,
                temperature: planet.atmosphere.temperatureProfile.troposphere,
                pressure: planet.atmosphere.surfacePressure,
                density: 0.19,      // kg/m³
                state: 'gas'
            }
        };
        
        // Heat sources
        this.heatSources = {
            // Primordial heat from formation
            primordialHeat: {
                power: 1.5e17, // W
                halfLife: 4.5e9 * 365.25 * 24 * 3600, // 4.5 billion years in seconds
                initialPower: 3e17
            },
            
            // Gravitational contraction (Kelvin-Helmholtz)
            gravitationalContraction: {
                power: 2.0e17, // W
                rate: -1e-8, // m/year contraction rate
                efficiency: 0.3
            },
            
            // Helium rain (differentiation)
            heliumRain: {
                power: 2.2e16, // W
                rate: 1e-12, // kg/s helium sinking rate
                releasedEnergy: 1.5e8 // J/kg
            },
            
            // Radioactive decay (minimal for gas giant)
            radioactiveDecay: {
                power: 1e15, // W (very small)
                isotopes: ['K40', 'U238', 'Th232']
            }
        };
        
        // Convection parameters
        this.convection = {
            rayleighNumber: 1e8, // Vigorous convection
            prandtlNumber: 0.8,
            nusseltNumber: 150,
            convectiveVelocity: 0.1, // m/s typical
            layerProperties: []
        };
        
        // Magnetic field generation (dynamo)
        this.dynamo = {
            active: true,
            region: 'metallicHydrogen', // Where metallic hydrogen exists
            electricalConductivity: 2e6, // S/m
            convectiveVelocity: 0.05, // m/s
            magneticReynolds: 1000,
            dipoleMoment: planet.magneticField.dipoleMoment
        };
        
        // Thermal evolution
        this.thermalHistory = {
            age: 4.6e9 * 365.25 * 24 * 3600, // 4.6 billion years in seconds
            coolingRate: -2e-13, // K/s
            luminosityEvolution: [],
            temperatureProfile: []
        };
        
        // Phase transitions
        this.phaseTransitions = {
            hydrogenMetallization: {
                pressure: 1.4e11, // Pa (140 GPa)
                temperature: 6000, // K
                depth: this.layers.metallicHydrogen.radius
            },
            
            heliumImmiscibility: {
                pressure: 1e10, // Pa (10 GPa)
                temperature: 1500, // K
                criticalPoint: true
            }
        };
        
        // Initialize core model
        this.initializeCore();
    }
    
    /**
     * Initialize the core thermal and structural model
     */
    initializeCore() {
        this.calculateLayerProperties();
        this.initializeConvection();
        this.calculateInitialTemperatureProfile();
        this.initializeMagneticDynamo();
    }
    
    /**
     * Calculate detailed properties for each layer
     */
    calculateLayerProperties() {
        Object.keys(this.layers).forEach(layerName => {
            const layer = this.layers[layerName];
            
            // Calculate mass of layer
            if (layerName === 'innerCore') {
                layer.volume = (4/3) * Math.PI * Math.pow(layer.radius * 1000, 3);
            } else {
                const outerRadius = layer.radius * 1000;
                const innerRadius = this.getPreviousLayerRadius(layerName) * 1000;
                layer.volume = (4/3) * Math.PI * (Math.pow(outerRadius, 3) - Math.pow(innerRadius, 3));
            }
            
            layer.mass = layer.density * layer.volume;
            
            // Calculate gravitational acceleration at layer boundary
            layer.gravity = this.calculateGravityAtRadius(layer.radius);
            
            // Calculate hydrostatic pressure
            layer.hydrostaticPressure = this.calculateHydrostaticPressure(layer.radius);
            
            // Calculate heat capacity
            layer.heatCapacity = this.calculateHeatCapacity(layer);
            
            // Calculate thermal conductivity
            layer.thermalConductivity = this.calculateThermalConductivity(layer);
        });
    }
    
    /**
     * Get radius of previous layer
     */
    getPreviousLayerRadius(layerName) {
        const layerOrder = ['innerCore', 'metallicHydrogen', 'molecularHydrogen', 'atmosphere'];
        const index = layerOrder.indexOf(layerName);
        if (index <= 0) return 0;
        return this.layers[layerOrder[index - 1]].radius;
    }
    
    /**
     * Calculate gravitational acceleration at given radius
     */
    calculateGravityAtRadius(radius) {
        let enclosedMass = 0;
        
        // Sum mass of all layers within radius
        Object.values(this.layers).forEach(layer => {
            if (layer.radius <= radius) {
                enclosedMass += layer.mass;
            } else {
                // Partial layer
                const fraction = Math.pow(radius / layer.radius, 3);
                enclosedMass += layer.mass * fraction;
            }
        });
        
        const G = 6.67430e-11; // m³/kg⋅s²
        const r = radius * 1000; // Convert to meters
        return G * enclosedMass / (r * r);
    }
    
    /**
     * Calculate hydrostatic pressure at given radius
     */
    calculateHydrostaticPressure(radius) {
        // Integrate pressure from surface inward
        let pressure = this.planet.atmosphere.surfacePressure;
        const steps = 100;
        const dr = (this.planet.radius - radius) / steps * 1000; // meters
        
        for (let i = 0; i < steps; i++) {
            const r = (this.planet.radius - i * dr / 1000) * 1000;
            const density = this.getDensityAtRadius(r / 1000);
            const gravity = this.calculateGravityAtRadius(r / 1000);
            
            pressure += density * gravity * dr;
        }
        
        return pressure;
    }
    
    /**
     * Get density at given radius
     */
    getDensityAtRadius(radius) {
        for (const layer of Object.values(this.layers)) {
            if (radius <= layer.radius) {
                return layer.density;
            }
        }
        return this.layers.atmosphere.density;
    }
    
    /**
     * Calculate heat capacity for layer
     */
    calculateHeatCapacity(layer) {
        // Simplified heat capacity based on composition
        let cp = 0;
        
        Object.entries(layer.composition).forEach(([component, fraction]) => {
            switch (component) {
                case 'H2':
                case 'metallicH2':
                    cp += fraction * 14300; // J/kg⋅K
                    break;
                case 'He':
                case 'metallicHe':
                    cp += fraction * 5200; // J/kg⋅K
                    break;
                case 'rock':
                    cp += fraction * 1000; // J/kg⋅K
                    break;
                case 'ice':
                    cp += fraction * 2100; // J/kg⋅K
                    break;
                default:
                    cp += fraction * 2000; // Default
            }
        });
        
        return cp;
    }
    
    /**
     * Calculate thermal conductivity for layer
     */
    calculateThermalConductivity(layer) {
        let k = 0;
        
        Object.entries(layer.composition).forEach(([component, fraction]) => {
            switch (component) {
                case 'metallicH2':
                case 'metallicHe':
                    k += fraction * 100; // W/m⋅K (high for metallic)
                    break;
                case 'H2':
                    k += fraction * 0.2; // W/m⋅K
                    break;
                case 'He':
                    k += fraction * 0.15; // W/m⋅K
                    break;
                case 'rock':
                    k += fraction * 3; // W/m⋅K
                    break;
                case 'ice':
                    k += fraction * 2.2; // W/m⋅K
                    break;
                default:
                    k += fraction * 0.5; // Default
            }
        });
        
        // Temperature and pressure corrections
        const T = layer.temperature;
        const P = layer.pressure;
        
        // Approximate temperature dependence
        k *= Math.pow(T / 300, 0.5);
        
        return k;
    }
    
    /**
     * Initialize convection model
     */
    initializeConvection() {
        Object.keys(this.layers).forEach(layerName => {
            const layer = this.layers[layerName];
            
            // Skip solid layers for convection
            if (layer.state === 'solid') return;
            
            // Calculate convective parameters
            const convectiveProperties = this.calculateConvectiveProperties(layer);
            
            this.convection.layerProperties.push({
                name: layerName,
                ...convectiveProperties
            });
        });
    }
    
    /**
     * Calculate convective properties for a layer
     */
    calculateConvectiveProperties(layer) {
        // Rayleigh number
        const g = layer.gravity;
        const alpha = 2e-4; // Thermal expansion coefficient (approximation)
        const deltaT = 1000; // Temperature difference across layer
        const d = layer.radius * 1000; // Layer thickness
        const nu = 1e-6; // Kinematic viscosity (approximation)
        const kappa = layer.thermalConductivity / (layer.density * layer.heatCapacity);
        
        const Ra = (g * alpha * deltaT * Math.pow(d, 3)) / (nu * kappa);
        
        // Nusselt number (heat transfer efficiency)
        const Nu = 0.27 * Math.pow(Ra, 0.25); // Simplified correlation
        
        // Convective velocity
        const v_conv = Math.sqrt(g * alpha * deltaT * d);
        
        return {
            rayleighNumber: Ra,
            nusseltNumber: Nu,
            convectiveVelocity: v_conv,
            thermalDiffusivity: kappa
        };
    }
    
    /**
     * Calculate initial temperature profile
     */
    calculateInitialTemperatureProfile() {
        const numPoints = 100;
        this.thermalHistory.temperatureProfile = [];
        
        for (let i = 0; i < numPoints; i++) {
            const radius = (i / (numPoints - 1)) * this.planet.radius;
            const temperature = this.calculateTemperatureAtRadius(radius);
            
            this.thermalHistory.temperatureProfile.push({
                radius: radius,
                temperature: temperature
            });
        }
    }
    
    /**
     * Calculate temperature at given radius
     */
    calculateTemperatureAtRadius(radius) {
        // Find which layer this radius is in
        let layer = this.layers.atmosphere;
        for (const [name, layerData] of Object.entries(this.layers)) {
            if (radius <= layerData.radius) {
                layer = layerData;
                break;
            }
        }
        
        // Interpolate temperature within layer
        const prevRadius = this.getPreviousLayerRadius(Object.keys(this.layers).find(key => this.layers[key] === layer));
        const fraction = (radius - prevRadius) / (layer.radius - prevRadius);
        
        // Temperature increases toward center (simplified)
        const surfaceTemp = this.planet.atmosphere.temperatureProfile.troposphere;
        const coreTemp = this.layers.innerCore.temperature;
        
        return surfaceTemp + (coreTemp - surfaceTemp) * Math.pow(1 - radius / this.planet.radius, 0.5);
    }
    
    /**
     * Initialize magnetic dynamo model
     */
    initializeMagneticDynamo() {
        const metalLayer = this.layers.metallicHydrogen;
        
        // Magnetic Reynolds number
        const velocity = metalLayer.convectiveVelocity || 0.05; // m/s
        const length = metalLayer.radius * 1000; // m
        const magneticDiffusivity = 1 / (4e-7 * Math.PI * metalLayer.electricalConductivity || 2e6);
        
        this.dynamo.magneticReynolds = velocity * length / magneticDiffusivity;
        
        // Dynamo is active if Rm > critical value (~10)
        this.dynamo.active = this.dynamo.magneticReynolds > 10;
        
        // Calculate magnetic field strength
        if (this.dynamo.active) {
            this.calculateMagneticField();
        }
    }
    
    /**
     * Calculate magnetic field generation
     */
    calculateMagneticField() {
        const metalLayer = this.layers.metallicHydrogen;
        
        // Approximate dipole moment from convective power
        const convectivePower = this.calculateConvectivePower(metalLayer);
        const efficiency = 1e-6; // Magnetic field generation efficiency
        
        // Magnetic energy
        const magneticEnergy = convectivePower * efficiency;
        const volume = metalLayer.volume;
        
        // Magnetic field strength
        const B = Math.sqrt(2 * 4e-7 * Math.PI * magneticEnergy / volume);
        
        // Update planet's magnetic field
        this.planet.magneticField.strength = B;
        this.planet.magneticField.dipoleMoment = B * Math.pow(metalLayer.radius * 1000, 3);
    }
    
    /**
     * Calculate convective power in a layer
     */
    calculateConvectivePower(layer) {
        const heatFlux = layer.thermalConductivity * 1000 / (layer.radius * 1000); // W/m²
        const area = 4 * Math.PI * Math.pow(layer.radius * 1000, 2);
        return heatFlux * area;
    }
    
    /**
     * Update heat sources over time
     */
    updateHeatSources(deltaTime) {
        // Primordial heat decay
        const primordial = this.heatSources.primordialHeat;
        const decayConstant = Math.log(2) / primordial.halfLife;
        primordial.power *= Math.exp(-decayConstant * deltaTime);
        
        // Gravitational contraction (decreases over time)
        const contraction = this.heatSources.gravitationalContraction;
        contraction.power *= (1 - 1e-10 * deltaTime); // Very slow decrease
        
        // Helium rain (may increase/decrease based on cooling)
        const heliumRain = this.heatSources.heliumRain;
        const coolingRate = this.thermalHistory.coolingRate;
        heliumRain.rate *= (1 + coolingRate * deltaTime * 1e-8);
        heliumRain.power = heliumRain.rate * heliumRain.releasedEnergy;
    }
    
    /**
     * Calculate total luminosity
     */
    getTotalLuminosity() {
        return Object.values(this.heatSources).reduce((total, source) => total + source.power, 0);
    }
    
    /**
     * Get current internal heat flux
     */
    getInternalHeatFlux() {
        const surfaceArea = 4 * Math.PI * Math.pow(this.planet.radius * 1000, 2);
        return this.getTotalLuminosity() / surfaceArea; // W/m²
    }
    
    /**
     * Update thermal evolution
     */
    updateThermalEvolution(deltaTime) {
        // Update heat sources
        this.updateHeatSources(deltaTime);
        
        // Update cooling rate
        const totalHeatCapacity = Object.values(this.layers).reduce((total, layer) => {
            return total + layer.mass * layer.heatCapacity;
        }, 0);
        
        const luminosity = this.getTotalLuminosity();
        this.thermalHistory.coolingRate = -luminosity / totalHeatCapacity;
        
        // Update temperature profile
        this.updateTemperatureProfile(deltaTime);
        
        // Update magnetic dynamo
        this.updateMagneticDynamo(deltaTime);
    }
    
    /**
     * Update temperature profile
     */
    updateTemperatureProfile(deltaTime) {
        this.thermalHistory.temperatureProfile.forEach(point => {
            // Cool based on luminosity and heat capacity
            const coolingRate = this.thermalHistory.coolingRate;
            point.temperature += coolingRate * deltaTime;
            
            // Minimum temperature constraint
            point.temperature = Math.max(point.temperature, 50); // K
        });
        
        // Update layer temperatures
        Object.values(this.layers).forEach(layer => {
            layer.temperature += this.thermalHistory.coolingRate * deltaTime * 0.1; // Slower for layers
        });
    }
    
    /**
     * Update magnetic dynamo
     */
    updateMagneticDynamo(deltaTime) {
        if (this.dynamo.active) {
            // Recalculate convective properties
            const metalLayer = this.layers.metallicHydrogen;
            const convProps = this.calculateConvectiveProperties(metalLayer);
            
            // Update convective velocity
            metalLayer.convectiveVelocity = convProps.convectiveVelocity;
            
            // Recalculate magnetic field
            this.calculateMagneticField();
        }
    }
    
    /**
     * Check for phase transitions
     */
    checkPhaseTransitions() {
        const transitions = [];
        
        // Check hydrogen metallization
        Object.values(this.layers).forEach(layer => {
            if (layer.pressure > this.phaseTransitions.hydrogenMetallization.pressure &&
                layer.temperature > this.phaseTransitions.hydrogenMetallization.temperature) {
                
                if (layer.composition.H2 && layer.state !== 'metallic') {
                    transitions.push({
                        type: 'hydrogenMetallization',
                        layer: layer,
                        newState: 'metallic'
                    });
                }
            }
        });
        
        return transitions;
    }
    
    /**
     * Get core structure information
     */
    getCoreInfo() {
        return {
            layers: this.layers,
            heatSources: this.heatSources,
            totalLuminosity: this.getTotalLuminosity(),
            internalHeatFlux: this.getInternalHeatFlux(),
            magneticDynamo: this.dynamo,
            coolingRate: this.thermalHistory.coolingRate,
            age: this.thermalHistory.age
        };
    }
    
    /**
     * Get temperature at specific radius
     */
    getTemperature(radius) {
        return this.calculateTemperatureAtRadius(radius);
    }
    
    /**
     * Get pressure at specific radius
     */
    getPressure(radius) {
        return this.calculateHydrostaticPressure(radius);
    }
    
    /**
     * Get density at specific radius
     */
    getDensity(radius) {
        return this.getDensityAtRadius(radius);
    }
    
    /**
     * Get composition at specific radius
     */
    getComposition(radius) {
        for (const layer of Object.values(this.layers)) {
            if (radius <= layer.radius) {
                return layer.composition;
            }
        }
        return this.layers.atmosphere.composition;
    }
    
    /**
     * Update core physics
     */
    update(deltaTime) {
        this.updateThermalEvolution(deltaTime);
        
        // Check for phase transitions
        const transitions = this.checkPhaseTransitions();
        if (transitions.length > 0) {
            console.log('Phase transitions detected:', transitions);
        }
    }
    
    /**
     * Export core state for saving
     */
    exportState() {
        return {
            heatSources: this.heatSources,
            thermalHistory: this.thermalHistory,
            dynamo: this.dynamo
        };
    }
    
    /**
     * Import core state from save
     */
    importState(state) {
        this.heatSources = { ...this.heatSources, ...state.heatSources };
        this.thermalHistory = { ...this.thermalHistory, ...state.thermalHistory };
        this.dynamo = { ...this.dynamo, ...state.dynamo };
    }
}
