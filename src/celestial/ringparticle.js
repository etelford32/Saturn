/**
 * RingParticle.js - Individual ring particle implementation
 * Handles physics, collisions, material properties, and interactions for Saturn's ring particles
 */

import * as THREE from 'three';
import { PlanetaryBody } from './PlanetaryBody.js';

export class RingParticle extends PlanetaryBody {
    constructor(config = {}) {
        // Initialize with minimal PlanetaryBody config
        const particleConfig = {
            name: config.name || `Particle_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
            mass: config.mass || 1e6, // kg (default ~1 meter ice particle)
            radius: config.radius || 1.0, // meters
            position: config.position || new THREE.Vector3(0, 0, 0),
            velocity: config.velocity || new THREE.Vector3(0, 0, 0),
            ...config
        };
        
        super(particleConfig);
        
        // Particle-specific properties
        this.particleId = this.generateParticleId();
        this.ringType = config.ringType || 'A'; // Which ring this particle belongs to
        this.generation = config.generation || 0; // How many collisions this particle has been through
        this.age = config.age || 0; // Seconds since particle creation/last major event
        
        // Physical composition and material properties
        this.composition = {
            ice: 0.95,    // Water ice fraction
            rock: 0.04,   // Rocky material fraction
            organics: 0.01, // Organic compounds fraction
            ...config.composition
        };
        
        // Calculate derived properties
        this.bulkDensity = this.calculateBulkDensity();
        this.porosity = config.porosity || this.calculatePorosity();
        this.solidDensity = this.bulkDensity / (1 - this.porosity);
        
        // Thermal properties
        this.thermal = {
            temperature: config.temperature || 80, // K (typical ring temperature)
            albedo: this.calculateAlbedo(),
            emissivity: config.emissivity || 0.95,
            thermalInertia: this.calculateThermalInertia(),
            heatCapacity: this.calculateHeatCapacity(),
            thermalConductivity: this.calculateThermalConductivity()
        };
        
        // Optical properties
        this.optical = {
            geometricAlbedo: this.thermal.albedo,
            phaseFunction: 'Henyey-Greenstein',
            asymmetryParameter: 0.3,
            extinctionEfficiency: 2.0,
            scatteringEfficiency: 1.8,
            absorptionEfficiency: 0.2,
            color: config.color || this.calculateColor()
        };
        
        // Electromagnetic properties (for spoke formation)
        this.electromagnetic = {
            charge: 0, // Coulombs
            chargeToMassRatio: 0, // C/kg
            magneticSusceptibility: this.calculateMagneticSusceptibility(),
            electricalResistivity: this.calculateElectricalResistivity(),
            dielectricConstant: this.calculateDielectricConstant()
        };
        
        // Collision properties
        this.collision = {
            restitutionCoefficient: config.restitution || this.calculateRestitution(),
            frictionCoefficient: config.friction || this.calculateFriction(),
            cohesionStrength: config.cohesion || this.calculateCohesion(), // N/m²
            tensileStrength: this.calculateTensileStrength(), // N/m²
            compressiveStrength: this.calculateCompressiveStrength(), // N/m²
            elasticModulus: this.calculateElasticModulus(), // Pa
            poissonRatio: config.poissonRatio || 0.25
        };
        
        // Shape and rotation
        this.shape = {
            type: config.shapeType || 'sphere', // sphere, ellipsoid, irregular
            aspectRatio: config.aspectRatio || 1.0, // For ellipsoids
            roughness: config.roughness || 0.1, // Surface roughness parameter
            orientation: config.orientation || new THREE.Euler(0, 0, 0),
            angularVelocity: config.angularVelocity || new THREE.Vector3(0, 0, 0)
        };
        
        // State tracking
        this.state = {
            phase: 'stable', // stable, colliding, fragmenting, aggregating, sublimating
            lastCollisionTime: 0,
            collisionCount: 0,
            sublimationRate: 0, // kg/s
            chargeAccumulation: 0,
            stressState: new THREE.Vector3(0, 0, 0), // Internal stress
            isTracked: config.tracked || false // For detailed analysis
        };
        
        // Environmental history
        this.environment = {
            radiationDose: 0, // J/kg accumulated
            micrometeoriteHits: 0,
            thermalCycles: 0,
            maximumTemperature: this.thermal.temperature,
            minimumTemperature: this.thermal.temperature,
            averageStress: 0
        };
        
        // Interaction radii
        this.interactionRadii = {
            collision: this.radius,
            gravity: this.radius * 5, // Effective gravitational influence
            electromagnetic: this.radius * 2,
            wake: this.calculateWakeRadius()
        };
        
        // Performance and optimization
        this.optimization = {
            updateFrequency: 1, // Update every N frames
            framesSinceUpdate: 0,
            isActive: true,
            neighborCount: 0,
            binIndex: -1 // Index in spatial bin
        };
        
        // Initialize particle
        this.initializeParticle();
    }
    
    /**
     * Generate unique particle ID
     */
    generateParticleId() {
        return `RP_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
    }
    
    /**
     * Initialize particle properties
     */
    initializeParticle() {
        this.updateDerivedProperties();
        this.initializeRotation();
        this.calculateInitialState();
    }
    
    /**
     * Calculate bulk density from composition
     */
    calculateBulkDensity() {
        let density = 0;
        
        // Component densities (kg/m³)
        const densities = {
            ice: 917,      // Water ice at ~80K
            rock: 2700,    // Silicate rock
            organics: 1200 // Organic compounds
        };
        
        Object.entries(this.composition).forEach(([component, fraction]) => {
            if (densities[component]) {
                density += fraction * densities[component];
            }
        });
        
        return density;
    }
    
    /**
     * Calculate porosity based on particle size and formation history
     */
    calculatePorosity() {
        // Smaller particles tend to be less porous
        const sizeFactor = Math.exp(-this.radius / 10); // meters
        
        // Base porosity depends on composition
        let basePorosity = 0.1; // 10% for pure ice
        if (this.composition.rock > 0.1) {
            basePorosity += 0.3; // Rocky particles more porous
        }
        
        // Generation affects porosity (more collisions = more compaction)
        const compactionFactor = Math.exp(-this.generation * 0.1);
        
        return Math.min(0.7, basePorosity * compactionFactor * (1 + sizeFactor));
    }
    
    /**
     * Calculate albedo based on composition and temperature
     */
    calculateAlbedo() {
        let albedo = 0;
        
        // Component albedos
        const albedos = {
            ice: 0.6,     // Clean water ice
            rock: 0.1,    // Dark rocky material
            organics: 0.05 // Very dark organics
        };
        
        Object.entries(this.composition).forEach(([component, fraction]) => {
            if (albedos[component]) {
                albedo += fraction * albedos[component];
            }
        });
        
        // Temperature dependence (warmer ice is darker)
        const tempFactor = Math.exp(-(this.thermal.temperature - 80) / 20);
        albedo *= tempFactor;
        
        // Age dependence (radiation darkening)
        const ageFactor = Math.exp(-this.environment.radiationDose / 1e6);
        albedo *= ageFactor;
        
        return Math.max(0.02, Math.min(0.9, albedo));
    }
    
    /**
     * Calculate thermal inertia
     */
    calculateThermalInertia() {
        const density = this.bulkDensity;
        const heatCapacity = this.calculateHeatCapacity();
        const conductivity = this.calculateThermalConductivity();
        
        return Math.sqrt(density * heatCapacity * conductivity);
    }
    
    /**
     * Calculate heat capacity
     */
    calculateHeatCapacity() {
        let cp = 0;
        
        // Component heat capacities (J/kg⋅K)
        const heatCapacities = {
            ice: 2100,     // Water ice
            rock: 1000,    // Silicate rock
            organics: 1500 // Organic compounds
        };
        
        Object.entries(this.composition).forEach(([component, fraction]) => {
            if (heatCapacities[component]) {
                cp += fraction * heatCapacities[component];
            }
        });
        
        // Temperature dependence
        cp *= (1 + 0.01 * (this.thermal.temperature - 80));
        
        return cp;
    }
    
    /**
     * Calculate thermal conductivity
     */
    calculateThermalConductivity() {
        let k = 0;
        
        // Component thermal conductivities (W/m⋅K)
        const conductivities = {
            ice: 2.2,      // Water ice
            rock: 3.0,     // Silicate rock
            organics: 0.5  // Organic compounds
        };
        
        Object.entries(this.composition).forEach(([component, fraction]) => {
            if (conductivities[component]) {
                k += fraction * conductivities[component];
            }
        });
        
        // Porosity effect (reduces conductivity)
        k *= Math.pow(1 - this.porosity, 1.5);
        
        return k;
    }
    
    /**
     * Calculate particle color based on composition and processing
     */
    calculateColor() {
        // Base colors for different compositions
        let color = new THREE.Color(1, 1, 1); // Start with white
        
        // Ice contribution (white to blue-white)
        const iceColor = new THREE.Color(0.9, 0.95, 1.0);
        color.lerp(iceColor, this.composition.ice);
        
        // Rock contribution (gray to brown)
        const rockColor = new THREE.Color(0.4, 0.3, 0.2);
        color.lerp(rockColor, this.composition.rock);
        
        // Organic contribution (dark brown to black)
        const organicColor = new THREE.Color(0.1, 0.05, 0.02);
        color.lerp(organicColor, this.composition.organics);
        
        // Radiation processing (darkening)
        const darkeningFactor = Math.exp(-this.environment.radiationDose / 1e5);
        color.multiplyScalar(darkeningFactor);
        
        return color.getHex();
    }
    
    /**
     * Calculate magnetic susceptibility
     */
    calculateMagneticSusceptibility() {
        // Most ring materials are diamagnetic (slightly negative susceptibility)
        let susceptibility = 0;
        
        const susceptibilities = {
            ice: -9e-6,    // Diamagnetic
            rock: 1e-4,    // Weakly paramagnetic (depends on iron content)
            organics: -5e-6 // Diamagnetic
        };
        
        Object.entries(this.composition).forEach(([component, fraction]) => {
            if (susceptibilities[component]) {
                susceptibility += fraction * susceptibilities[component];
            }
        });
        
        return susceptibility;
    }
    
    /**
     * Calculate electrical resistivity
     */
    calculateElectricalResistivity() {
        // Pure ice is highly resistive, but impurities reduce resistivity
        let resistivity = 1e15; // Ohm⋅m for pure ice
        
        // Rock content reduces resistivity
        if (this.composition.rock > 0.01) {
            resistivity *= Math.exp(-this.composition.rock * 10);
        }
        
        // Temperature dependence
        resistivity *= Math.exp(5000 / this.thermal.temperature - 5000 / 273);
        
        return resistivity;
    }
    
    /**
     * Calculate dielectric constant
     */
    calculateDielectricConstant() {
        let dielectric = 1;
        
        const dielectrics = {
            ice: 3.2,      // Water ice
            rock: 8.0,     // Silicate rock
            organics: 2.5  // Organic compounds
        };
        
        Object.entries(this.composition).forEach(([component, fraction]) => {
            if (dielectrics[component]) {
                dielectric += fraction * (dielectrics[component] - 1);
            }
        });
        
        // Porosity effect
        dielectric = 1 + (dielectric - 1) * (1 - this.porosity);
        
        return dielectric;
    }
    
    /**
     * Calculate restitution coefficient
     */
    calculateRestitution() {
        // Ice particles are somewhat elastic
        let restitution = 0.3;
        
        // Rock content reduces elasticity
        restitution *= Math.exp(-this.composition.rock * 2);
        
        // Temperature dependence (colder = more brittle)
        restitution *= Math.min(1, this.thermal.temperature / 100);
        
        // Size dependence (smaller particles more elastic)
        restitution *= Math.exp(-this.radius / 5);
        
        return Math.max(0.05, Math.min(0.8, restitution));
    }
    
    /**
     * Calculate friction coefficient
     */
    calculateFriction() {
        let friction = 0.3; // Base friction for ice
        
        // Rock content increases friction
        friction += this.composition.rock * 0.5;
        
        // Temperature dependence
        if (this.thermal.temperature > 130) {
            // Surface melting increases friction
            friction *= 1.5;
        }
        
        // Surface roughness effect
        friction *= (1 + this.shape.roughness);
        
        return Math.max(0.1, Math.min(1.0, friction));
    }
    
    /**
     * Calculate cohesion strength
     */
    calculateCohesion() {
        // Van der Waals forces and sintering
        let cohesion = 100; // Pa base cohesion
        
        // Ice content increases cohesion (sintering)
        cohesion *= (1 + this.composition.ice * 5);
        
        // Temperature dependence (sintering is stronger at higher temps)
        if (this.thermal.temperature > 100) {
            cohesion *= Math.pow(this.thermal.temperature / 100, 2);
        }
        
        // Size dependence (surface area effect)
        cohesion *= Math.exp(-this.radius / 10);
        
        return cohesion;
    }
    
    /**
     * Calculate tensile strength
     */
    calculateTensileStrength() {
        // Tensile strength much lower than compressive for granular materials
        const compressive = this.calculateCompressiveStrength();
        return compressive * 0.1; // Rule of thumb for ice/rock
    }
    
    /**
     * Calculate compressive strength
     */
    calculateCompressiveStrength() {
        let strength = 5e6; // Pa for solid ice
        
        // Porosity reduces strength
        strength *= Math.pow(1 - this.porosity, 2);
        
        // Rock content can increase strength
        strength *= (1 + this.composition.rock * 2);
        
        // Temperature dependence
        strength *= Math.exp(-(this.thermal.temperature - 80) / 50);
        
        return strength;
    }
    
    /**
     * Calculate elastic modulus
     */
    calculateElasticModulus() {
        let modulus = 9e9; // Pa for solid ice
        
        // Porosity effect
        modulus *= Math.pow(1 - this.porosity, 2);
        
        // Composition effect
        modulus *= (1 + this.composition.rock * 5); // Rock is stiffer
        
        // Temperature dependence
        modulus *= Math.exp(-(this.thermal.temperature - 80) / 100);
        
        return modulus;
    }
    
    /**
     * Calculate wake radius for gravitational wakes
     */
    calculateWakeRadius() {
        // Toomre Q parameter for stability
        const Q = 2.0; // Stable disk value
        const sigma = 10; // kg/m² surface density (typical)
        const kappa = 2 * Math.PI / (this.getOrbitalPeriod() || 86400); // Epicyclic frequency
        
        // Wake wavelength
        const lambda = 2 * Math.PI * Q * sigma / (this.bulkDensity * kappa);
        
        return lambda / (2 * Math.PI); // Wake radius
    }
    
    /**
     * Initialize particle rotation
     */
    initializeRotation() {
        // Random initial orientation
        this.shape.orientation = new THREE.Euler(
            Math.random() * 2 * Math.PI,
            Math.random() * 2 * Math.PI,
            Math.random() * 2 * Math.PI
        );
        
        // Random angular velocity (slow tumbling)
        const maxSpin = 0.1; // rad/s
        this.shape.angularVelocity = new THREE.Vector3(
            (Math.random() - 0.5) * maxSpin,
            (Math.random() - 0.5) * maxSpin,
            (Math.random() - 0.5) * maxSpin
        );
    }
    
    /**
     * Calculate initial state
     */
    calculateInitialState() {
        this.state.phase = 'stable';
        this.state.lastCollisionTime = 0;
        this.state.collisionCount = 0;
        this.electromagnetic.charge = 0;
        this.electromagnetic.chargeToMassRatio = 0;
    }
    
    /**
     * Update derived properties
     */
    updateDerivedProperties() {
        this.bulkDensity = this.calculateBulkDensity();
        this.thermal.albedo = this.calculateAlbedo();
        this.thermal.thermalInertia = this.calculateThermalInertia();
        this.optical.color = this.calculateColor();
        this.interactionRadii.wake = this.calculateWakeRadius();
    }
    
    /**
     * Apply thermal evolution
     */
    updateThermal(deltaTime, solarFlux = 0, infraredFlux = 0) {
        // Solar heating
        const solarHeating = solarFlux * (1 - this.thermal.albedo) * Math.PI * this.radius * this.radius;
        
        // Infrared heating from Saturn
        const irHeating = infraredFlux * this.thermal.emissivity * Math.PI * this.radius * this.radius;
        
        // Radiative cooling (Stefan-Boltzmann)
        const stefanBoltzmann = 5.67e-8; // W/m²⋅K⁴
        const radiativeCooling = 4 * Math.PI * this.radius * this.radius * this.thermal.emissivity * 
                               stefanBoltzmann * Math.pow(this.thermal.temperature, 4);
        
        // Net heat flow
        const netHeat = solarHeating + irHeating - radiativeCooling;
        
        // Temperature change
        const mass = this.mass;
        const deltaT = netHeat * deltaTime / (mass * this.thermal.heatCapacity);
        
        this.thermal.temperature += deltaT;
        
        // Update temperature extremes
        this.environment.maximumTemperature = Math.max(this.environment.maximumTemperature, this.thermal.temperature);
        this.environment.minimumTemperature = Math.min(this.environment.minimumTemperature, this.thermal.temperature);
        
        // Check for sublimation (ice)
        if (this.thermal.temperature > 150 && this.composition.ice > 0) {
            this.updateSublimation(deltaTime);
        }
    }
    
    /**
     * Update sublimation process
     */
    updateSublimation(deltaTime) {
        // Simplified sublimation rate (Hertz-Knudsen equation)
        const vaporPressure = this.getVaporPressure(this.thermal.temperature);
        const molecularMass = 18e-3; // kg/mol for water
        const gasConstant = 8.314; // J/mol⋅K
        const surfaceArea = 4 * Math.PI * this.radius * this.radius;
        
        const sublimationRate = vaporPressure * surfaceArea * 
                              Math.sqrt(molecularMass / (2 * Math.PI * gasConstant * this.thermal.temperature));
        
        this.state.sublimationRate = sublimationRate;
        
        // Mass loss
        const massLoss = sublimationRate * deltaTime;
        if (massLoss > 0 && this.mass > massLoss) {
            this.mass -= massLoss;
            this.composition.ice = Math.max(0, this.composition.ice - massLoss / this.mass);
            
            // Recalculate radius
            this.radius = Math.cbrt(3 * this.mass / (4 * Math.PI * this.bulkDensity));
            
            this.updateDerivedProperties();
        }
    }
    
    /**
     * Get vapor pressure of ice
     */
    getVaporPressure(temperature) {
        // Antoine equation for ice
        const A = 9.550426;
        const B = 5723.265;
        const C = 3.53068;
        
        const logP = A - B / (temperature + C);
        return Math.pow(10, logP); // Pa
    }
    
    /**
     * Update electromagnetic state
     */
    updateElectromagnetic(deltaTime, plasmaEnvironment = {}) {
        // Charge accumulation from plasma
        const electronFlux = plasmaEnvironment.electronFlux || 0;
        const ionFlux = plasmaEnvironment.ionFlux || 0;
        
        // Surface area
        const area = 4 * Math.PI * this.radius * this.radius;
        
        // Charge accumulation rate
        const chargeRate = (ionFlux - electronFlux) * area * 1.6e-19; // Elementary charge
        
        // Apply charging
        this.electromagnetic.charge += chargeRate * deltaTime;
        this.electromagnetic.chargeToMassRatio = this.electromagnetic.charge / this.mass;
        
        // Charge dissipation through conductivity
        const dissipationRate = this.electromagnetic.charge / 
                              (this.electromagnetic.electricalResistivity * this.electromagnetic.dielectricConstant);
        
        this.electromagnetic.charge *= Math.exp(-dissipationRate * deltaTime);
        
        // Spoke formation (electrostatic levitation)
        if (Math.abs(this.electromagnetic.chargeToMassRatio) > 1e-6 && this.radius < 0.1) {
            this.state.phase = 'electrostatic_levitation';
        }
    }
    
    /**
     * Update rotation state
     */
    updateRotation(deltaTime) {
        // Apply angular velocity to orientation
        this.shape.orientation.x += this.shape.angularVelocity.x * deltaTime;
        this.shape.orientation.y += this.shape.angularVelocity.y * deltaTime;
        this.shape.orientation.z += this.shape.angularVelocity.z * deltaTime;
        
        // Damping due to tidal effects (very small)
        const dampingFactor = 1 - 1e-10 * deltaTime;
        this.shape.angularVelocity.multiplyScalar(dampingFactor);
    }
    
    /**
     * Handle collision with another particle
     */
    processCollision(otherParticle, impactVelocity, contactPoint) {
        this.state.collisionCount++;
        this.state.lastCollisionTime = Date.now();
        this.generation = Math.max(this.generation, otherParticle.generation) + 1;
        
        const collisionEnergy = 0.5 * this.getReducedMass(otherParticle) * impactVelocity * impactVelocity;
        
        // Determine collision outcome
        if (collisionEnergy > this.getFragmentationThreshold(otherParticle)) {
            this.state.phase = 'fragmenting';
            return this.fragmentationOutcome(otherParticle, collisionEnergy);
        } else if (impactVelocity < this.getAggregationThreshold(otherParticle)) {
            this.state.phase = 'aggregating';
            return this.aggregationOutcome(otherParticle);
        } else {
            this.state.phase = 'bouncing';
            return this.bounceOutcome(otherParticle, impactVelocity);
        }
    }
    
    /**
     * Get reduced mass for collision calculations
     */
    getReducedMass(otherParticle) {
        return (this.mass * otherParticle.mass) / (this.mass + otherParticle.mass);
    }
    
    /**
     * Get fragmentation threshold
     */
    getFragmentationThreshold(otherParticle) {
        // Based on tensile strength and particle size
        const characteristicSize = Math.min(this.radius, otherParticle.radius);
        const characteristicStrength = Math.min(this.collision.tensileStrength, otherParticle.collision.tensileStrength);
        
        return characteristicStrength * characteristicSize * characteristicSize * characteristicSize;
    }
    
    /**
     * Get aggregation threshold
     */
    getAggregationThreshold(otherParticle) {
        // Based on cohesion strength
        const combinedCohesion = Math.sqrt(this.collision.cohesionStrength * otherParticle.collision.cohesionStrength);
        const contactArea = Math.PI * Math.pow(Math.min(this.radius, otherParticle.radius), 2);
        
        return Math.sqrt(2 * combinedCohesion * contactArea / this.getReducedMass(otherParticle));
    }
    
    /**
     * Calculate fragmentation outcome
     */
    fragmentationOutcome(otherParticle, energy) {
        // Simple fragmentation model
        const totalMass = this.mass + otherParticle.mass;
        const fragmentCount = Math.floor(2 + Math.sqrt(energy / 1e6));
        
        return {
            type: 'fragmentation',
            fragmentCount: fragmentCount,
            totalMass: totalMass,
            energyDissipated: energy * 0.9 // Most energy goes into creating new surfaces
        };
    }
    
    /**
     * Calculate aggregation outcome
     */
    aggregationOutcome(otherParticle) {
        const totalMass = this.mass + otherParticle.mass;
        const combinedRadius = Math.cbrt(Math.pow(this.radius, 3) + Math.pow(otherParticle.radius, 3));
        
        // Weighted average composition
        const newComposition = {};
        Object.keys(this.composition).forEach(component => {
            const weighted = (this.composition[component] * this.mass + 
                            (otherParticle.composition[component] || 0) * otherParticle.mass) / totalMass;
            newComposition[component] = weighted;
        });
        
        return {
            type: 'aggregation',
            newMass: totalMass,
            newRadius: combinedRadius,
            newComposition: newComposition
        };
    }
    
    /**
     * Calculate bounce outcome
     */
    bounceOutcome(otherParticle, impactVelocity) {
        return {
            type: 'bounce',
            energyLoss: impactVelocity * impactVelocity * (1 - this.collision.restitutionCoefficient) / 2
        };
    }
    
    /**
     * Update particle age and evolution
     */
    updateEvolution(deltaTime) {
        this.age += deltaTime;
        
        // Radiation damage accumulation
        const radiationRate = 1e-10; // J/kg⋅s (very rough estimate)
        this.environment.radiationDose += radiationRate * deltaTime;
        
        // Micrometeorite bombardment (statistical)
        const impactProbability = 1e-12 * deltaTime * 4 * Math.PI * this.radius * this.radius;
        if (Math.random() < impactProbability) {
            this.environment.micrometeoriteHits++;
            this.processMinorImpact();
        }
        
        // Thermal cycling
        const thermalCycleRate = 1 / (this.getOrbitalPeriod() || 86400);
        this.environment.thermalCycles += thermalCycleRate * deltaTime;
        
        // Update properties based on aging
        if (this.environment.radiationDose > 1e6) {
            this.thermal.albedo *= 0.99999; // Gradual darkening
        }
    }
    
    /**
     * Process minor impact event
     */
    processMinorImpact() {
        // Small mass loss
        const massLoss = this.mass * 1e-8;
        this.mass = Math.max(this.mass - massLoss, this.mass * 0.01);
        
        // Slight heating
        this.thermal.temperature += 1;
        
        // Surface roughening
        this.shape.roughness = Math.min(1, this.shape.roughness * 1.01);
    }
    
    /**
     * Check if particle is stable
     */
    isStable() {
        return this.state.phase === 'stable' && 
               this.mass > 0 && 
               this.radius > this.physics?.minimumParticleSize;
    }
    
    /**
     * Get particle information
     */
    getInfo() {
        return {
            id: this.particleId,
            ringType: this.ringType,
            mass: this.mass,
            radius: this.radius,
            composition: this.composition,
            temperature: this.thermal.temperature,
            albedo: this.thermal.albedo,
            charge: this.electromagnetic.charge,
            age: this.age,
            generation: this.generation,
            collisions: this.state.collisionCount,
            phase: this.state.phase,
            position: this.position.toArray(),
            velocity: this.velocity.toArray()
        };
    }
    
    /**
     * Main particle update
     */
    update(deltaTime, environment = {}) {
        if (!this.isStable()) return;
        
        this.optimization.framesSinceUpdate++;
        
        // Adaptive update frequency based on activity
        if (this.optimization.framesSinceUpdate >= this.optimization.updateFrequency) {
            // Update thermal state
            this.updateThermal(deltaTime, environment.solarFlux, environment.infraredFlux);
            
            // Update electromagnetic state
            this.updateElectromagnetic(deltaTime, environment.plasma);
            
            // Update rotation
            this.updateRotation(deltaTime);
            
            // Update evolution
            this.updateEvolution(deltaTime);
            
            // Reset frame counter
            this.optimization.framesSinceUpdate = 0;
        }
        
        // Always update position/velocity (handled by parent class)
        // this.updatePhysics(deltaTime) would be called by ring system
    }
    
    /**
     * Export particle state
     */
    exportState() {
        return {
            ...super.exportState(),
            particleId: this.particleId,
            ringType: this.ringType,
            generation: this.generation,
            age: this.age,
            composition: this.composition,
            thermal: this.thermal,
            electromagnetic: this.electromagnetic,
            collision: this.collision,
            shape: {
                ...this.shape,
                orientation: this.shape.orientation.toArray(),
                angularVelocity: this.shape.angularVelocity.toArray()
            },
            state: this.state,
            environment: this.environment
        };
    }
    
    /**
     * Import particle state
     */
    importState(state) {
        super.importState(state);
        
        if (state.particleId) this.particleId = state.particleId;
        if (state.ringType) this.ringType = state.ringType;
        if (state.generation !== undefined) this.generation = state.generation;
        if (state.age !== undefined) this.age = state.age;
        if (state.composition) this.composition = { ...state.composition };
        if (state.thermal) this.thermal = { ...this.thermal, ...state.thermal };
        if (state.electromagnetic) this.electromagnetic = { ...this.electromagnetic, ...state.electromagnetic };
        if (state.collision) this.collision = { ...this.collision, ...state.collision };
        if (state.shape) {
            this.shape = { ...this.shape, ...state.shape };
            if (state.shape.orientation) this.shape.orientation.fromArray(state.shape.orientation);
            if (state.shape.angularVelocity) this.shape.angularVelocity.fromArray(state.shape.angularVelocity);
        }
        if (state.state) this.state = { ...this.state, ...state.state };
        if (state.environment) this.environment = { ...this.environment, ...state.environment };
        
        this.updateDerivedProperties();
    }
}
