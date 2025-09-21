/**
 * PlanetaryBody.js - Base class for all celestial objects in the Saturn system
 * Provides fundamental properties and methods for planets, moons, and ring particles
 */

import * as THREE from 'three';

export class PlanetaryBody {
    constructor(config = {}) {
        // Physical properties
        this.name = config.name || 'Unnamed Body';
        this.mass = config.mass || 1.0; // kg
        this.radius = config.radius || 1.0; // km (will be scaled for rendering)
        this.density = config.density || this.calculateDensity(); // kg/m³
        this.gravParameter = config.gravParameter || this.calculateGravParameter(); // GM (m³/s²)
        
        // Orbital elements (Keplerian)
        this.semiMajorAxis = config.semiMajorAxis || 0; // km
        this.eccentricity = config.eccentricity || 0;
        this.inclination = config.inclination || 0; // radians
        this.longitudeOfAscendingNode = config.longitudeOfAscendingNode || 0; // radians
        this.argumentOfPeriapsis = config.argumentOfPeriapsis || 0; // radians
        this.meanAnomalyAtEpoch = config.meanAnomalyAtEpoch || 0; // radians
        this.epochTime = config.epochTime || 0; // seconds since J2000
        
        // Current state vectors
        this.position = new THREE.Vector3(0, 0, 0); // km
        this.velocity = new THREE.Vector3(0, 0, 0); // km/s
        this.acceleration = new THREE.Vector3(0, 0, 0); // km/s²
        
        // Rotational properties
        this.rotationPeriod = config.rotationPeriod || 86400; // seconds (1 day default)
        this.axialTilt = config.axialTilt || 0; // radians
        this.rotationPhase = config.rotationPhase || 0; // current rotation angle
        
        // Physical constants
        this.albedo = config.albedo || 0.1; // geometric albedo
        this.temperature = config.temperature || 0; // Kelvin
        this.atmosphericPressure = config.atmosphericPressure || 0; // Pascals
        
        // Rendering properties
        this.color = config.color || 0xFFFFFF;
        this.renderScale = config.renderScale || 1.0;
        this.visible = config.visible !== undefined ? config.visible : true;
        
        // Three.js mesh (will be created by subclasses)
        this.mesh = null;
        this.orbitLine = null;
        
        // Hierarchy relationships
        this.parent = null;
        this.children = [];
        
        // Physics tracking
        this.forces = [];
        this.lastUpdateTime = 0;
        
        // Scientific metadata
        this.metadata = {
            discoveryDate: config.discoveryDate || null,
            discoveredBy: config.discoveredBy || null,
            classification: config.classification || 'Unknown',
            composition: config.composition || {},
            magneticField: config.magneticField || null
        };
    }
    
    /**
     * Calculate density if not provided
     */
    calculateDensity() {
        if (this.mass && this.radius) {
            const volumeM3 = (4/3) * Math.PI * Math.pow(this.radius * 1000, 3);
            return this.mass / volumeM3;
        }
        return 1000; // Default density (water)
    }
    
    /**
     * Calculate gravitational parameter (GM)
     */
    calculateGravParameter() {
        const G = 6.67430e-11; // m³/kg⋅s²
        return G * this.mass;
    }
    
    /**
     * Get current mean motion (radians per second)
     */
    getMeanMotion() {
        if (this.semiMajorAxis <= 0) return 0;
        const a_m = this.semiMajorAxis * 1000; // Convert to meters
        const parentMu = this.parent ? this.parent.gravParameter : 3.7931187e16; // Saturn's GM
        return Math.sqrt(parentMu / Math.pow(a_m, 3));
    }
    
    /**
     * Convert orbital elements to state vectors (position and velocity)
     */
    orbitalElementsToStateVectors(time) {
        if (!this.parent || this.semiMajorAxis <= 0) {
            return { position: this.position.clone(), velocity: this.velocity.clone() };
        }
        
        const n = this.getMeanMotion();
        const M = this.meanAnomalyAtEpoch + n * (time - this.epochTime);
        
        // Solve Kepler's equation for eccentric anomaly
        const E = this.solveKeplersEquation(M, this.eccentricity);
        
        // True anomaly
        const nu = 2 * Math.atan2(
            Math.sqrt(1 + this.eccentricity) * Math.sin(E/2),
            Math.sqrt(1 - this.eccentricity) * Math.cos(E/2)
        );
        
        // Distance from focus
        const r = this.semiMajorAxis * (1 - this.eccentricity * Math.cos(E));
        
        // Position in orbital plane
        const cos_nu = Math.cos(nu);
        const sin_nu = Math.sin(nu);
        const x_orb = r * cos_nu;
        const y_orb = r * sin_nu;
        
        // Velocity in orbital plane
        const sqrt_mu_a = Math.sqrt(this.parent.gravParameter / (this.semiMajorAxis * 1000));
        const sqrt_1_e2 = Math.sqrt(1 - this.eccentricity * this.eccentricity);
        const vx_orb = -sqrt_mu_a * sin_nu / sqrt_1_e2;
        const vy_orb = sqrt_mu_a * (this.eccentricity + cos_nu) / sqrt_1_e2;
        
        // Rotation matrices for orbital plane to inertial frame
        const cos_W = Math.cos(this.longitudeOfAscendingNode);
        const sin_W = Math.sin(this.longitudeOfAscendingNode);
        const cos_w = Math.cos(this.argumentOfPeriapsis);
        const sin_w = Math.sin(this.argumentOfPeriapsis);
        const cos_i = Math.cos(this.inclination);
        const sin_i = Math.sin(this.inclination);
        
        // Transform to inertial frame
        const position = new THREE.Vector3(
            (cos_W * cos_w - sin_W * sin_w * cos_i) * x_orb + (-cos_W * sin_w - sin_W * cos_w * cos_i) * y_orb,
            (sin_W * cos_w + cos_W * sin_w * cos_i) * x_orb + (-sin_W * sin_w + cos_W * cos_w * cos_i) * y_orb,
            (sin_w * sin_i) * x_orb + (cos_w * sin_i) * y_orb
        );
        
        const velocity = new THREE.Vector3(
            (cos_W * cos_w - sin_W * sin_w * cos_i) * vx_orb + (-cos_W * sin_w - sin_W * cos_w * cos_i) * vy_orb,
            (sin_W * cos_w + cos_W * sin_w * cos_i) * vx_orb + (-sin_W * sin_w + cos_W * cos_w * cos_i) * vy_orb,
            (sin_w * sin_i) * vx_orb + (cos_w * sin_i) * vy_orb
        ).multiplyScalar(1e-3); // Convert m/s to km/s
        
        return { position, velocity };
    }
    
    /**
     * Solve Kepler's equation using Newton-Raphson method
     */
    solveKeplersEquation(M, e, tolerance = 1e-12, maxIterations = 100) {
        let E = M; // Initial guess
        
        for (let i = 0; i < maxIterations; i++) {
            const f = E - e * Math.sin(E) - M;
            const fp = 1 - e * Math.cos(E);
            const delta = f / fp;
            E -= delta;
            
            if (Math.abs(delta) < tolerance) {
                break;
            }
        }
        
        return E;
    }
    
    /**
     * Add a force to be applied during next physics update
     */
    addForce(force) {
        this.forces.push(force);
    }
    
    /**
     * Clear all accumulated forces
     */
    clearForces() {
        this.forces.length = 0;
        this.acceleration.set(0, 0, 0);
    }
    
    /**
     * Calculate gravitational acceleration from another body
     */
    calculateGravitationalAcceleration(otherBody) {
        const separation = otherBody.position.clone().sub(this.position);
        const distance = separation.length() * 1000; // Convert km to m
        
        if (distance === 0) return new THREE.Vector3(0, 0, 0);
        
        const acceleration = separation.normalize().multiplyScalar(
            otherBody.gravParameter / (distance * distance) / 1000 // Convert back to km/s²
        );
        
        return acceleration;
    }
    
    /**
     * Update physics state using numerical integration
     */
    updatePhysics(deltaTime, integrator = 'leapfrog') {
        this.lastUpdateTime += deltaTime;
        
        switch (integrator) {
            case 'leapfrog':
                this.leapfrogIntegration(deltaTime);
                break;
            case 'rk4':
                this.rungeKutta4Integration(deltaTime);
                break;
            default:
                this.eulerIntegration(deltaTime);
        }
        
        this.clearForces();
    }
    
    /**
     * Leapfrog integration (symplectic, energy-conserving)
     */
    leapfrogIntegration(dt) {
        // Calculate total acceleration from all forces
        const totalAcceleration = new THREE.Vector3(0, 0, 0);
        this.forces.forEach(force => totalAcceleration.add(force));
        
        // Update velocity (kick)
        this.velocity.add(totalAcceleration.clone().multiplyScalar(dt * 0.5));
        
        // Update position (drift)
        this.position.add(this.velocity.clone().multiplyScalar(dt));
        
        // Update velocity (kick)
        this.velocity.add(totalAcceleration.clone().multiplyScalar(dt * 0.5));
        
        this.acceleration = totalAcceleration;
    }
    
    /**
     * Runge-Kutta 4th order integration (higher accuracy)
     */
    rungeKutta4Integration(dt) {
        const pos0 = this.position.clone();
        const vel0 = this.velocity.clone();
        
        const totalAcceleration = new THREE.Vector3(0, 0, 0);
        this.forces.forEach(force => totalAcceleration.add(force));
        
        // k1
        const k1v = totalAcceleration.clone().multiplyScalar(dt);
        const k1r = vel0.clone().multiplyScalar(dt);
        
        // k2
        const k2v = totalAcceleration.clone().multiplyScalar(dt);
        const k2r = vel0.clone().add(k1v.clone().multiplyScalar(0.5)).multiplyScalar(dt);
        
        // k3
        const k3v = totalAcceleration.clone().multiplyScalar(dt);
        const k3r = vel0.clone().add(k2v.clone().multiplyScalar(0.5)).multiplyScalar(dt);
        
        // k4
        const k4v = totalAcceleration.clone().multiplyScalar(dt);
        const k4r = vel0.clone().add(k3v).multiplyScalar(dt);
        
        // Final update
        this.velocity = vel0.add(
            k1v.add(k2v.multiplyScalar(2)).add(k3v.multiplyScalar(2)).add(k4v).multiplyScalar(1/6)
        );
        this.position = pos0.add(
            k1r.add(k2r.multiplyScalar(2)).add(k3r.multiplyScalar(2)).add(k4r).multiplyScalar(1/6)
        );
        
        this.acceleration = totalAcceleration;
    }
    
    /**
     * Simple Euler integration (less accurate but faster)
     */
    eulerIntegration(dt) {
        const totalAcceleration = new THREE.Vector3(0, 0, 0);
        this.forces.forEach(force => totalAcceleration.add(force));
        
        this.velocity.add(totalAcceleration.clone().multiplyScalar(dt));
        this.position.add(this.velocity.clone().multiplyScalar(dt));
        
        this.acceleration = totalAcceleration;
    }
    
    /**
     * Update rotational state
     */
    updateRotation(deltaTime) {
        const angularVelocity = (2 * Math.PI) / this.rotationPeriod;
        this.rotationPhase += angularVelocity * deltaTime;
        this.rotationPhase %= (2 * Math.PI);
        
        if (this.mesh) {
            this.mesh.rotation.y = this.rotationPhase;
            this.mesh.rotation.z = this.axialTilt;
        }
    }
    
    /**
     * Update render position based on physics position
     */
    updateRenderPosition(scaleFactor = 1) {
        if (this.mesh) {
            // Apply scale factor for visualization (distances in space are huge)
            this.mesh.position.copy(this.position.clone().multiplyScalar(scaleFactor));
        }
    }
    
    /**
     * Calculate orbital period in seconds
     */
    getOrbitalPeriod() {
        if (!this.parent || this.semiMajorAxis <= 0) return 0;
        const a_m = this.semiMajorAxis * 1000;
        return 2 * Math.PI * Math.sqrt(Math.pow(a_m, 3) / this.parent.gravParameter);
    }
    
    /**
     * Calculate escape velocity from surface
     */
    getEscapeVelocity() {
        return Math.sqrt(2 * this.gravParameter / (this.radius * 1000)) / 1000; // km/s
    }
    
    /**
     * Calculate surface gravity
     */
    getSurfaceGravity() {
        return this.gravParameter / Math.pow(this.radius * 1000, 2); // m/s²
    }
    
    /**
     * Get Hill sphere radius (sphere of gravitational influence)
     */
    getHillSphereRadius() {
        if (!this.parent) return Infinity;
        const massRatio = this.mass / this.parent.mass;
        return this.semiMajorAxis * Math.pow(massRatio / 3, 1/3);
    }
    
    /**
     * Check if point is within Roche limit
     */
    isWithinRocheLimit(distance) {
        if (!this.parent) return false;
        const densityRatio = this.parent.density / this.density;
        const rocheLimit = 2.44 * this.parent.radius * Math.pow(densityRatio, 1/3);
        return distance < rocheLimit;
    }
    
    /**
     * Calculate tidal acceleration at given position
     */
    calculateTidalAcceleration(position) {
        if (!this.parent) return new THREE.Vector3(0, 0, 0);
        
        const toCenter = this.parent.position.clone().sub(position);
        const toCenterMag = toCenter.length() * 1000; // Convert to m
        const toParent = this.parent.position.clone().sub(this.position);
        const toParentMag = toParent.length() * 1000;
        
        // Tidal acceleration = GM * (r_vector/|r|³ - R_vector/|R|³)
        const tidalAccel = toCenter.normalize().multiplyScalar(this.parent.gravParameter / Math.pow(toCenterMag, 3))
            .sub(toParent.normalize().multiplyScalar(this.parent.gravParameter / Math.pow(toParentMag, 3)));
        
        return tidalAccel.multiplyScalar(1e-3); // Convert to km/s²
    }
    
    /**
     * Export state for serialization
     */
    exportState() {
        return {
            name: this.name,
            position: this.position.toArray(),
            velocity: this.velocity.toArray(),
            rotationPhase: this.rotationPhase,
            lastUpdateTime: this.lastUpdateTime
        };
    }
    
    /**
     * Import state from serialization
     */
    importState(state) {
        this.position.fromArray(state.position);
        this.velocity.fromArray(state.velocity);
        this.rotationPhase = state.rotationPhase;
        this.lastUpdateTime = state.lastUpdateTime;
    }
    
    /**
     * Dispose of Three.js resources
     */
    dispose() {
        if (this.mesh) {
            if (this.mesh.geometry) this.mesh.geometry.dispose();
            if (this.mesh.material) {
                if (Array.isArray(this.mesh.material)) {
                    this.mesh.material.forEach(material => material.dispose());
                } else {
                    this.mesh.material.dispose();
                }
            }
        }
        
        if (this.orbitLine) {
            if (this.orbitLine.geometry) this.orbitLine.geometry.dispose();
            if (this.orbitLine.material) this.orbitLine.material.dispose();
        }
    }
}
