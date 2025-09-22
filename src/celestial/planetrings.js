/**
 * PlanetRings.js - Saturn's ring system implementation
 * Manages the complex ring structure, particle dynamics, and shepherd moon interactions
 */

import * as THREE from 'three';
import { RingParticle } from './RingParticle.js';

export class PlanetRings {
    constructor(planet, config = {}) {
        this.planet = planet;
        
        // Ring system configuration
        this.config = {
            particleCount: config.particleCount || 50000,
            maxParticles: config.maxParticles || 200000,
            enableCollisions: config.enableCollisions !== false,
            enableSelfGravity: config.enableSelfGravity !== false,
            enableShepherdMoons: config.enableShepherdMoons !== false,
            qualityLevel: config.qualityLevel || 'medium', // low, medium, high, ultra
            renderMode: config.renderMode || 'particles', // particles, instanced, compute
            ...config
        };
        
        // Saturn's ring structure (accurate distances in Saturn radii)
        this.ringStructure = {
            D: {
                name: 'D Ring',
                innerRadius: 1.11, // 66,900 km
                outerRadius: 1.235, // 74,500 km
                opticalDepth: 1e-6,
                particleDensity: 1e-8, // kg/m³
                color: 0xD3D3D3,
                composition: { ice: 0.95, rock: 0.05 }
            },
            C: {
                name: 'C Ring (Crepe Ring)',
                innerRadius: 1.235, // 74,500 km
                outerRadius: 1.525, // 92,000 km
                opticalDepth: 0.05,
                particleDensity: 1e-6,
                color: 0xC0C0C0,
                composition: { ice: 0.90, rock: 0.10 },
                features: ['Maxwell Gap', 'Bond Gap']
            },
            B: {
                name: 'B Ring',
                innerRadius: 1.525, // 92,000 km
                outerRadius: 1.950, // 117,580 km
                opticalDepth: 2.5,
                particleDensity: 5e-5,
                color: 0xE0E0E0,
                composition: { ice: 0.95, rock: 0.05 },
                features: ['Spokes'] // Electrostatic phenomena
            },
            CassiniDivision: {
                name: 'Cassini Division',
                innerRadius: 1.950, // 117,580 km
                outerRadius: 2.025, // 122,170 km
                opticalDepth: 0.01,
                particleDensity: 1e-7,
                color: 0x808080,
                composition: { ice: 0.85, rock: 0.15 },
                resonance: '2:1 Mimas'
            },
            A: {
                name: 'A Ring',
                innerRadius: 2.025, // 122,170 km
                outerRadius: 2.267, // 136,780 km
                opticalDepth: 0.4,
                particleDensity: 1e-5,
                color: 0xF0F0F0,
                composition: { ice: 0.93, rock: 0.07 },
                features: ['Encke Gap', 'Keeler Gap']
            },
            F: {
                name: 'F Ring',
                innerRadius: 2.326, // 140,180 km
                outerRadius: 2.329, // 140,180 km (very narrow)
                opticalDepth: 0.1,
                particleDensity: 1e-6,
                color: 0xFFFFFF,
                composition: { ice: 0.98, rock: 0.02 },
                shepherds: ['Prometheus', 'Pandora']
            },
            E: {
                name: 'E Ring',
                innerRadius: 3.0, // 181,000 km
                outerRadius: 8.0, // 483,000 km
                opticalDepth: 1e-7,
                particleDensity: 1e-9,
                color: 0xE6F3FF,
                composition: { ice: 0.99, organics: 0.01 },
                source: 'Enceladus geysers'
            }
        };
        
        // Particle system
        this.particles = [];
        this.particlePool = [];
        this.spatialHash = new Map(); // For collision detection
        this.hashGridSize = 2000; // meters
        
        // Ring bins for efficient simulation
        this.ringBins = [];
        this.binWidth = 1000; // meters
        this.minRadius = this.ringStructure.D.innerRadius * this.planet.radius * 1000;
        this.maxRadius = this.ringStructure.E.outerRadius * this.planet.radius * 1000;
        
        // Shepherd moons
        this.shepherdMoons = [];
        this.moonInteractionRadius = 50000; // meters
        
        // Physics parameters
        this.physics = {
            gravitationalParameter: this.planet.gravParameter,
            selfGravityStrength: 6.67430e-11, // G
            dampingCoefficient: 0.1, // Energy loss in collisions
            minimumParticleSize: 0.01, // meters
            maximumParticleSize: 10.0, // meters
            sizeDistributionIndex: -3.0, // Power law index
            rotationRate: 2 * Math.PI / this.planet.rotationPeriod
        };
        
        // Collision parameters
        this.collision = {
            restitutionCoefficient: 0.1, // Inelastic collisions
            frictionCoefficient: 0.3,
            cohesionStrength: 1.0, // N/m² (weak)
            fragmentationThreshold: 100.0, // m/s collision speed
            aggregationProbability: 0.01
        };
        
        // Optical properties
        this.optics = {
            particleAlbedo: 0.6,
            phaseFunction: 'Henyey-Greenstein',
            asymmetryParameter: 0.3,
            extinctionEfficiency: 2.0,
            scatteringEfficiency: 1.8
        };
        
        // Rendering components
        this.mesh = null;
        this.instancedMesh = null;
        this.computeShader = null;
        this.material = null;
        
        // Performance monitoring
        this.performance = {
            collisionChecks: 0,
            collisionEvents: 0,
            gravityCalculations: 0,
            frameTime: 0
        };
        
        // Toomre Q parameter for self-gravity stability
        this.toomreQ = 2.0; // Stable value
        
        this.initializeRings();
    }
    
    /**
     * Initialize the ring system
     */
    initializeRings() {
        this.createRingBins();
        this.generateParticles();
        this.setupSpatialHashing();
        this.createRenderingComponents();
        this.initializeShepherdMoons();
    }
    
    /**
     * Create ring bins for efficient simulation
     */
    createRingBins() {
        const numBins = Math.ceil((this.maxRadius - this.minRadius) / this.binWidth);
        
        for (let i = 0; i < numBins; i++) {
            const innerR = this.minRadius + i * this.binWidth;
            const outerR = Math.min(innerR + this.binWidth, this.maxRadius);
            const centerR = (innerR + outerR) / 2;
            
            // Determine which ring this bin belongs to
            const ringType = this.getRingTypeAtRadius(centerR / (this.planet.radius * 1000));
            const ringData = this.ringStructure[ringType];
            
            this.ringBins.push({
                index: i,
                innerRadius: innerR,
                outerRadius: outerR,
                centerRadius: centerR,
                particles: [],
                ringType: ringType,
                opticalDepth: ringData.opticalDepth,
                surfaceDensity: this.calculateSurfaceDensity(centerR, ringData),
                orbitalPeriod: this.calculateOrbitalPeriod(centerR),
                meanMotion: Math.sqrt(this.physics.gravitationalParameter / Math.pow(centerR, 3)),
                epicyclicFrequency: 0, // Will be calculated
                verticalFrequency: 0
            });
        }
        
        // Calculate frequencies for each bin
        this.calculateBinFrequencies();
    }
    
    /**
     * Get ring type at given radius (in Saturn radii)
     */
    getRingTypeAtRadius(radius) {
        for (const [type, ring] of Object.entries(this.ringStructure)) {
            if (radius >= ring.innerRadius && radius <= ring.outerRadius) {
                return type;
            }
        }
        return 'E'; // Default to E ring for outer regions
    }
    
    /**
     * Calculate surface density for a ring
     */
    calculateSurfaceDensity(radius, ringData) {
        // Use optical depth to estimate surface density
        const particleRadius = 1.0; // meters (average)
        const crossSection = Math.PI * particleRadius * particleRadius;
        const particleDensity = 917; // kg/m³ (ice)
        
        // τ = σ * Q / (ρ * r) where σ is surface density
        return ringData.opticalDepth * particleDensity * particleRadius / this.optics.extinctionEfficiency;
    }
    
    /**
     * Calculate orbital period at radius
     */
    calculateOrbitalPeriod(radius) {
        return 2 * Math.PI * Math.sqrt(Math.pow(radius, 3) / this.physics.gravitationalParameter);
    }
    
    /**
     * Calculate epicyclic and vertical frequencies
     */
    calculateBinFrequencies() {
        this.ringBins.forEach(bin => {
            const r = bin.centerRadius;
            const n = bin.meanMotion;
            
            // Epicyclic frequency (radial oscillations)
            bin.epicyclicFrequency = n * Math.sqrt(1 + (3/4) * Math.pow(this.planet.radius * 1000 / r, 2) * this.planet.J2);
            
            // Vertical frequency (out-of-plane oscillations)
            bin.verticalFrequency = n * Math.sqrt(1 - (3/4) * Math.pow(this.planet.radius * 1000 / r, 2) * this.planet.J2);
        });
    }
    
    /**
     * Generate ring particles
     */
    generateParticles() {
        const particlesPerRing = this.config.particleCount / Object.keys(this.ringStructure).length;
        
        Object.entries(this.ringStructure).forEach(([ringType, ringData]) => {
            if (ringType === 'E') {
                // E ring gets fewer particles due to its large size and low density
                this.generateRingParticles(ringType, ringData, particlesPerRing * 0.1);
            } else {
                this.generateRingParticles(ringType, ringData, particlesPerRing);
            }
        });
        
        console.log(`Generated ${this.particles.length} ring particles`);
    }
    
    /**
     * Generate particles for a specific ring
     */
    generateRingParticles(ringType, ringData, count) {
        const planetRadius = this.planet.radius * 1000; // meters
        const innerRadius = ringData.innerRadius * planetRadius;
        const outerRadius = ringData.outerRadius * planetRadius;
        
        for (let i = 0; i < count; i++) {
            // Random radius (weighted by area)
            const u = Math.random();
            const radius = Math.sqrt(u * (outerRadius * outerRadius - innerRadius * innerRadius) + innerRadius * innerRadius);
            
            // Random azimuth
            const azimuth = Math.random() * 2 * Math.PI;
            
            // Particle size from power-law distribution
            const minSize = this.physics.minimumParticleSize;
            const maxSize = this.physics.maximumParticleSize;
            const sizeRandom = Math.random();
            const size = minSize * Math.pow(maxSize / minSize, sizeRandom);
            
            // Position
            const position = new THREE.Vector3(
                radius * Math.cos(azimuth),
                (Math.random() - 0.5) * size * 10, // Small vertical dispersion
                radius * Math.sin(azimuth)
            );
            
            // Keplerian velocity with perturbations
            const orbitalSpeed = Math.sqrt(this.physics.gravitationalParameter / radius);
            const velocity = new THREE.Vector3(
                -orbitalSpeed * Math.sin(azimuth) + (Math.random() - 0.5) * 10,
                (Math.random() - 0.5) * 1, // Small vertical velocity
                orbitalSpeed * Math.cos(azimuth) + (Math.random() - 0.5) * 10
            );
            
            // Create particle
            const particle = new RingParticle({
                position: position,
                velocity: velocity,
                radius: size,
                mass: this.calculateParticleMass(size, ringData.composition),
                composition: ringData.composition,
                ringType: ringType,
                color: ringData.color
            });
            
            this.particles.push(particle);
            
            // Add to appropriate bin
            const binIndex = this.getBinIndex(radius);
            if (binIndex >= 0 && binIndex < this.ringBins.length) {
                this.ringBins[binIndex].particles.push(particle);
                particle.binIndex = binIndex;
            }
        }
    }
    
    /**
     * Calculate particle mass from size and composition
     */
    calculateParticleMass(radius, composition) {
        const volume = (4/3) * Math.PI * Math.pow(radius, 3);
        let density = 0;
        
        // Weighted average density
        Object.entries(composition).forEach(([material, fraction]) => {
            switch (material) {
                case 'ice':
                    density += fraction * 917; // kg/m³
                    break;
                case 'rock':
                    density += fraction * 2700; // kg/m³
                    break;
                case 'organics':
                    density += fraction * 1200; // kg/m³
                    break;
            }
        });
        
        return volume * density;
    }
    
    /**
     * Get bin index for a given radius
     */
    getBinIndex(radius) {
        return Math.floor((radius - this.minRadius) / this.binWidth);
    }
    
    /**
     * Setup spatial hashing for collision detection
     */
    setupSpatialHashing() {
        this.spatialHash.clear();
        
        this.particles.forEach(particle => {
            const hash = this.getParticleHash(particle);
            if (!this.spatialHash.has(hash)) {
                this.spatialHash.set(hash, []);
            }
            this.spatialHash.get(hash).push(particle);
        });
    }
    
    /**
     * Get spatial hash for particle
     */
    getParticleHash(particle) {
        const gridX = Math.floor(particle.position.x / this.hashGridSize);
        const gridY = Math.floor(particle.position.y / this.hashGridSize);
        const gridZ = Math.floor(particle.position.z / this.hashGridSize);
        return `${gridX},${gridY},${gridZ}`;
    }
    
    /**
     * Initialize shepherd moons
     */
    initializeShepherdMoons() {
        // F-ring shepherds
        this.shepherdMoons.push({
            name: 'Prometheus',
            radius: 139380000, // meters from Saturn center
            mass: 1.595e17, // kg
            influence: 5000, // meters
            type: 'inner_shepherd'
        });
        
        this.shepherdMoons.push({
            name: 'Pandora',
            radius: 141700000, // meters
            mass: 1.371e17, // kg
            influence: 5000, // meters
            type: 'outer_shepherd'
        });
        
        // A-ring shepherds
        this.shepherdMoons.push({
            name: 'Pan',
            radius: 133583000, // meters (Encke Gap)
            mass: 4.95e15, // kg
            influence: 2000, // meters
            type: 'gap_shepherd'
        });
        
        this.shepherdMoons.push({
            name: 'Daphnis',
            radius: 136505000, // meters (Keeler Gap)
            mass: 7.7e13, // kg
            influence: 1000, // meters
            type: 'gap_shepherd'
        });
    }
    
    /**
     * Create rendering components
     */
    createRenderingComponents() {
        switch (this.config.renderMode) {
            case 'instanced':
                this.createInstancedMesh();
                break;
            case 'compute':
                this.createComputeShader();
                break;
            default:
                this.createParticleSystem();
        }
    }
    
    /**
     * Create standard particle system
     */
    createParticleSystem() {
        const geometry = new THREE.BufferGeometry();
        const positions = new Float32Array(this.particles.length * 3);
        const colors = new Float32Array(this.particles.length * 3);
        const sizes = new Float32Array(this.particles.length);
        
        this.particles.forEach((particle, i) => {
            positions[i * 3] = particle.position.x;
            positions[i * 3 + 1] = particle.position.y;
            positions[i * 3 + 2] = particle.position.z;
            
            const color = new THREE.Color(particle.color);
            colors[i * 3] = color.r;
            colors[i * 3 + 1] = color.g;
            colors[i * 3 + 2] = color.b;
            
            sizes[i] = particle.radius * 2; // Diameter for rendering
        });
        
        geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        geometry.setAttribute('size', new THREE.BufferAttribute(sizes, 1));
        
        const material = new THREE.PointsMaterial({
            size: 10,
            vertexColors: true,
            transparent: true,
            opacity: 0.8,
            blending: THREE.AdditiveBlending,
            sizeAttenuation: true
        });
        
        this.mesh = new THREE.Points(geometry, material);
        this.material = material;
    }
    
    /**
     * Create instanced mesh for better performance
     */
    createInstancedMesh() {
        const geometry = new THREE.SphereGeometry(1, 8, 6);
        const material = new THREE.MeshBasicMaterial({ 
            transparent: true, 
            opacity: 0.6,
            vertexColors: true
        });
        
        this.instancedMesh = new THREE.InstancedMesh(geometry, material, this.particles.length);
        
        const matrix = new THREE.Matrix4();
        const color = new THREE.Color();
        
        this.particles.forEach((particle, i) => {
            matrix.makeTranslation(particle.position.x, particle.position.y, particle.position.z);
            matrix.multiplyScalar(particle.radius);
            this.instancedMesh.setMatrixAt(i, matrix);
            
            color.setHex(particle.color);
            this.instancedMesh.setColorAt(i, color);
        });
        
        this.instancedMesh.instanceMatrix.needsUpdate = true;
        this.instancedMesh.instanceColor.needsUpdate = true;
        
        this.mesh = this.instancedMesh;
    }
    
    /**
     * Update particle physics
     */
    updatePhysics(deltaTime) {
        const startTime = performance.now();
        
        // Clear performance counters
        this.performance.collisionChecks = 0;
        this.performance.collisionEvents = 0;
        this.performance.gravityCalculations = 0;
        
        // Update particle orbits
        this.updateOrbitalMotion(deltaTime);
        
        // Apply shepherd moon perturbations
        if (this.config.enableShepherdMoons) {
            this.applyShepherdMoonEffects(deltaTime);
        }
        
        // Handle collisions
        if (this.config.enableCollisions) {
            this.handleCollisions(deltaTime);
        }
        
        // Apply self-gravity
        if (this.config.enableSelfGravity) {
            this.applySelfGravity(deltaTime);
        }
        
        // Update spatial hashing
        this.updateSpatialHashing();
        
        this.performance.frameTime = performance.now() - startTime;
    }
    
    /**
     * Update orbital motion of particles
     */
    updateOrbitalMotion(deltaTime) {
        this.particles.forEach(particle => {
            // Gravitational acceleration from Saturn
            const r = particle.position.length();
            if (r > 0) {
                const accel = -this.physics.gravitationalParameter / (r * r * r);
                const acceleration = particle.position.clone().multiplyScalar(accel);
                
                // Add J2 oblateness effect
                if (this.planet.J2) {
                    const j2Accel = this.calculateJ2Acceleration(particle.position);
                    acceleration.add(j2Accel);
                }
                
                // Update velocity and position
                particle.velocity.add(acceleration.clone().multiplyScalar(deltaTime));
                particle.position.add(particle.velocity.clone().multiplyScalar(deltaTime));
                
                // Update bin if particle moved
                const newBinIndex = this.getBinIndex(r);
                if (newBinIndex !== particle.binIndex && newBinIndex >= 0 && newBinIndex < this.ringBins.length) {
                    this.moveParticleToBin(particle, newBinIndex);
                }
            }
        });
    }
    
    /**
     * Calculate J2 oblateness acceleration
     */
    calculateJ2Acceleration(position) {
        const r = position.length();
        const Re = this.planet.radius * 1000; // meters
        const J2 = this.planet.J2 || 1.6298e-2;
        const mu = this.physics.gravitationalParameter;
        
        const x = position.x;
        const y = position.y;
        const z = position.z;
        const r2 = r * r;
        const r5 = r2 * r2 * r;
        const Re2 = Re * Re;
        
        const factor = -1.5 * J2 * mu * Re2 / r5;
        
        const ax = factor * x * (1 - 5 * z * z / r2);
        const ay = factor * y * (1 - 5 * z * z / r2);
        const az = factor * z * (3 - 5 * z * z / r2);
        
        return new THREE.Vector3(ax, ay, az);
    }
    
    /**
     * Move particle to new bin
     */
    moveParticleToBin(particle, newBinIndex) {
        // Remove from old bin
        if (particle.binIndex >= 0 && particle.binIndex < this.ringBins.length) {
            const oldBin = this.ringBins[particle.binIndex];
            const index = oldBin.particles.indexOf(particle);
            if (index > -1) {
                oldBin.particles.splice(index, 1);
            }
        }
        
        // Add to new bin
        this.ringBins[newBinIndex].particles.push(particle);
        particle.binIndex = newBinIndex;
    }
    
    /**
     * Apply shepherd moon gravitational effects
     */
    applyShepherdMoonEffects(deltaTime) {
        this.shepherdMoons.forEach(moon => {
            // Calculate moon position (simplified circular orbit)
            const moonAngle = (Date.now() / 1000) * Math.sqrt(this.physics.gravitationalParameter / Math.pow(moon.radius, 3));
            const moonPos = new THREE.Vector3(
                moon.radius * Math.cos(moonAngle),
                0,
                moon.radius * Math.sin(moonAngle)
            );
            
            // Find particles within influence
            this.particles.forEach(particle => {
                const distance = particle.position.distanceTo(moonPos);
                
                if (distance < moon.influence) {
                    // Calculate gravitational perturbation
                    const direction = moonPos.clone().sub(particle.position).normalize();
                    const force = 6.67430e-11 * moon.mass / (distance * distance);
                    const acceleration = direction.multiplyScalar(force / particle.mass);
                    
                    particle.velocity.add(acceleration.multiplyScalar(deltaTime));
                    
                    // Gap maintenance for gap shepherds
                    if (moon.type === 'gap_shepherd') {
                        this.maintainGap(particle, moon, distance);
                    }
                }
            });
        });
    }
    
    /**
     * Maintain gaps created by shepherd moons
     */
    maintainGap(particle, moon, distance) {
        const gapRadius = moon.radius;
        const particleRadius = particle.position.length();
        
        // Repel particles from gap center
        if (Math.abs(particleRadius - gapRadius) < moon.influence) {
            const repulsionStrength = 1e-6; // m/s²
            const direction = particleRadius > gapRadius ? 1 : -1;
            
            particle.velocity.add(
                particle.position.clone().normalize().multiplyScalar(direction * repulsionStrength)
            );
        }
    }
    
    /**
     * Handle particle collisions
     */
    handleCollisions(deltaTime) {
        // Use spatial hashing for efficient collision detection
        for (const [hash, particles] of this.spatialHash) {
            if (particles.length > 1) {
                for (let i = 0; i < particles.length; i++) {
                    for (let j = i + 1; j < particles.length; j++) {
                        this.checkCollision(particles[i], particles[j], deltaTime);
                        this.performance.collisionChecks++;
                    }
                }
            }
        }
    }
    
    /**
     * Check collision between two particles
     */
    checkCollision(p1, p2, deltaTime) {
        const distance = p1.position.distanceTo(p2.position);
        const minDistance = p1.radius + p2.radius;
        
        if (distance < minDistance) {
            this.resolveCollision(p1, p2);
            this.performance.collisionEvents++;
        }
    }
    
    /**
     * Resolve collision between particles
     */
    resolveCollision(p1, p2) {
        // Collision normal
        const normal = p1.position.clone().sub(p2.position).normalize();
        
        // Relative velocity
        const relativeVelocity = p1.velocity.clone().sub(p2.velocity);
        const velocityAlongNormal = relativeVelocity.dot(normal);
        
        // Do not resolve if velocities are separating
        if (velocityAlongNormal > 0) return;
        
        // Collision impulse
        const e = this.collision.restitutionCoefficient;
        const impulse = -(1 + e) * velocityAlongNormal / (1/p1.mass + 1/p2.mass);
        
        // Apply impulse
        const impulseVector = normal.clone().multiplyScalar(impulse);
        p1.velocity.add(impulseVector.clone().multiplyScalar(1/p1.mass));
        p2.velocity.sub(impulseVector.clone().multiplyScalar(1/p2.mass));
        
        // Separate particles
        const overlap = (p1.radius + p2.radius) - p1.position.distanceTo(p2.position);
        const separation = normal.clone().multiplyScalar(overlap / 2);
        p1.position.add(separation);
        p2.position.sub(separation);
        
        // Check for fragmentation or aggregation
        const collisionSpeed = Math.abs(velocityAlongNormal);
        if (collisionSpeed > this.collision.fragmentationThreshold) {
            this.fragmentParticles(p1, p2, collisionSpeed);
        } else if (collisionSpeed < 1 && Math.random() < this.collision.aggregationProbability) {
            this.aggregateParticles(p1, p2);
        }
    }
    
    /**
     * Fragment particles in high-energy collisions
     */
    fragmentParticles(p1, p2, collisionSpeed) {
        // Simple fragmentation: create smaller particles
        const fragments = [];
        const totalMass = p1.mass + p2.mass;
        const numFragments = Math.floor(2 + Math.random() * 4);
        
        for (let i = 0; i < numFragments; i++) {
            const fragmentMass = totalMass / numFragments;
            const fragmentRadius = Math.cbrt(fragmentMass / (p1.mass / Math.pow(p1.radius, 3)));
            
            if (fragmentRadius > this.physics.minimumParticleSize) {
                const fragment = new RingParticle({
                    position: p1.position.clone().add(
                        new THREE.Vector3(
                            (Math.random() - 0.5) * p1.radius,
                            (Math.random() - 0.5) * p1.radius,
                            (Math.random() - 0.5) * p1.radius
                        )
                    ),
                    velocity: p1.velocity.clone().add(
                        new THREE.Vector3(
                            (Math.random() - 0.5) * collisionSpeed * 0.1,
                            (Math.random() - 0.5) * collisionSpeed * 0.1,
                            (Math.random() - 0.5) * collisionSpeed * 0.1
                        )
                    ),
                    radius: fragmentRadius,
                    mass: fragmentMass,
                    composition: p1.composition,
                    ringType: p1.ringType,
                    color: p1.color
                });
                
                fragments.push(fragment);
            }
        }
        
        // Remove original particles and add fragments
        if (fragments.length > 0) {
            this.removeParticle(p1);
            this.removeParticle(p2);
            fragments.forEach(fragment => this.addParticle(fragment));
        }
    }
    
    /**
     * Aggregate particles in low-energy collisions
     */
    aggregateParticles(p1, p2) {
        // Create new particle from combination
        const totalMass = p1.mass + p2.mass;
        const newRadius = Math.cbrt(totalMass / (p1.mass / Math.pow(p1.radius, 3)));
        
        if (newRadius < this.physics.maximumParticleSize) {
            const centerOfMass = p1.position.clone().multiplyScalar(p1.mass)
                .add(p2.position.clone().multiplyScalar(p2.mass))
                .divideScalar(totalMass);
            
            const combinedVelocity = p1.velocity.clone().multiplyScalar(p1.mass)
                .add(p2.velocity.clone().multiplyScalar(p2.mass))
                .divideScalar(totalMass);
            
            const aggregate = new RingParticle({
                position: centerOfMass,
                velocity: combinedVelocity,
                radius: newRadius,
                mass: totalMass,
                composition: p1.composition, // Assume same composition
                ringType: p1.ringType,
                color: p1.color
            });
            
            this.removeParticle(p1);
            this.removeParticle(p2);
            this.addParticle(aggregate);
        }
    }
    
    /**
     * Apply self-gravity between particles
     */
    applySelfGravity(deltaTime) {
        // Only apply to dense regions to save computation
        this.ringBins.forEach(bin => {
            if (bin.particles.length > 5) {
                this.calculateBinSelfGravity(bin, deltaTime);
            }
        });
    }
    
    /**
     * Calculate self-gravity within a bin
     */
    calculateBinSelfGravity(bin, deltaTime) {
        const particles = bin.particles;
        
        for (let i = 0; i < particles.length; i++) {
            let totalForce = new THREE.Vector3(0, 0, 0);
            
            for (let j = 0; j < particles.length; j++) {
                if (i !== j) {
                    const force = this.calculateGravitationalForce(particles[i], particles[j]);
                    totalForce.add(force);
                    this.performance.gravityCalculations++;
                }
            }
            
            // Apply force
            const acceleration = totalForce.divideScalar(particles[i].mass);
            particles[i].velocity.add(acceleration.multiplyScalar(deltaTime));
        }
    }
    
    /**
     * Calculate gravitational force between two particles
     */
    calculateGravitationalForce(p1, p2) {
        const separation = p2.position.clone().sub(p1.position);
        const distance = separation.length();
        
        if (distance === 0) return new THREE.Vector3(0, 0, 0);
        
        const forceMagnitude = this.physics.selfGravityStrength * p1.mass * p2.mass / (distance * distance);
        return separation.normalize().multiplyScalar(forceMagnitude);
    }
    
    /**
     * Update spatial hashing
     */
    updateSpatialHashing() {
        this.spatialHash.clear();
        
        this.particles.forEach(particle => {
            const hash = this.getParticleHash(particle);
            if (!this.spatialHash.has(hash)) {
                this.spatialHash.set(hash, []);
            }
            this.spatialHash.get(hash).push(particle);
        });
    }
    
    /**
     * Add particle to system
     */
    addParticle(particle) {
        if (this.particles.length < this.config.maxParticles) {
            this.particles.push(particle);
            
            const binIndex = this.getBinIndex(particle.position.length());
            if (binIndex >= 0 && binIndex < this.ringBins.length) {
                this.ringBins[binIndex].particles.push(particle);
                particle.binIndex = binIndex;
            }
        }
    }
    
    /**
     * Remove particle from system
     */
    removeParticle(particle) {
        const index = this.particles.indexOf(particle);
        if (index > -1) {
            this.particles.splice(index, 1);
            
            // Remove from bin
            if (particle.binIndex >= 0 && particle.binIndex < this.ringBins.length) {
                const bin = this.ringBins[particle.binIndex];
                const binIndex = bin.particles.indexOf(particle);
                if (binIndex > -1) {
                    bin.particles.splice(binIndex, 1);
                }
            }
        }
    }
    
    /**
     * Calculate optical depth at radius
     */
    getOpticalDepth(radius) {
        const binIndex = this.getBinIndex(radius);
        if (binIndex >= 0 && binIndex < this.ringBins.length) {
            return this.ringBins[binIndex].opticalDepth;
        }
        return 0;
    }
    
    /**
     * Get ring information
     */
    getRingInfo() {
        return {
            totalParticles: this.particles.length,
            activeBins: this.ringBins.filter(bin => bin.particles.length > 0).length,
            totalBins: this.ringBins.length,
            performance: this.performance,
            ringStructure: this.ringStructure,
            shepherdMoons: this.shepherdMoons
        };
    }
    
    /**
     * Update ring system
     */
    update(deltaTime, timeScale = 1) {
        const scaledDelta = deltaTime * timeScale;
        
        // Update physics
        this.updatePhysics(scaledDelta);
        
        // Update rendering
        this.updateRendering();
    }
    
    /**
     * Update rendering components
     */
    updateRendering() {
        if (this.mesh && this.mesh.geometry && this.mesh.geometry.attributes.position) {
            const positions = this.mesh.geometry.attributes.position.array;
            const colors = this.mesh.geometry.attributes.color.array;
            
            this.particles.forEach((particle, i) => {
                if (i * 3 + 2 < positions.length) {
                    positions[i * 3] = particle.position.x;
                    positions[i * 3 + 1] = particle.position.y;
                    positions[i * 3 + 2] = particle.position.z;
                }
            });
            
            this.mesh.geometry.attributes.position.needsUpdate = true;
        }
        
        // Update instanced mesh if using instancing
        if (this.instancedMesh) {
            const matrix = new THREE.Matrix4();
            
            this.particles.forEach((particle, i) => {
                matrix.makeTranslation(particle.position.x, particle.position.y, particle.position.z);
                matrix.multiplyScalar(particle.radius);
                this.instancedMesh.setMatrixAt(i, matrix);
            });
            
            this.instancedMesh.instanceMatrix.needsUpdate = true;
        }
    }
    
    /**
     * Add to scene
     */
    addToScene(scene) {
        if (this.mesh) {
            scene.add(this.mesh);
        }
    }
    
    /**
     * Remove from scene
     */
    removeFromScene(scene) {
        if (this.mesh) {
            scene.remove(this.mesh);
        }
    }
    
    /**
     * Export ring state
     */
    exportState() {
        return {
            particles: this.particles.map(p => p.exportState()),
            performance: this.performance,
            config: this.config
        };
    }
    
    /**
     * Import ring state
     */
    importState(state) {
        // This would restore particle positions and velocities
        console.log('Importing ring state...');
    }
    
    /**
     * Dispose of resources
     */
    dispose() {
        if (this.mesh) {
            if (this.mesh.geometry) this.mesh.geometry.dispose();
            if (this.mesh.material) this.mesh.material.dispose();
        }
        
        if (this.instancedMesh) {
            this.instancedMesh.geometry.dispose();
            this.instancedMesh.material.dispose();
        }
        
        this.particles.length = 0;
        this.spatialHash.clear();
    }
}
