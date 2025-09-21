/**
 * Planet.js - Saturn-specific implementation
 * Handles Saturn's unique properties, atmosphere, magnetosphere, and complex geometry
 */

import * as THREE from 'three';
import { PlanetaryBody } from './PlanetaryBody.js';

export class Planet extends PlanetaryBody {
    constructor(config = {}) {
        // Saturn-specific default configuration
        const saturnConfig = {
            name: 'Saturn',
            mass: 5.6834e26, // kg
            radius: 58232, // km (equatorial radius)
            polarRadius: 54364, // km (polar radius due to oblateness)
            density: 687, // kg/m³ (less dense than water!)
            rotationPeriod: 38520, // seconds (10.7 hours)
            axialTilt: THREE.MathUtils.degToRad(26.73), // radians
            albedo: 0.47,
            temperature: 95, // K (-178°C)
            gravParameter: 3.7931187e16, // m³/s²
            color: 0xFAD5A5,
            ...config
        };
        
        super(saturnConfig);
        
        // Oblateness parameters (gravitational harmonics)
        this.J2 = 1.6298e-2; // Second zonal harmonic
        this.J4 = -9.35e-4; // Fourth zonal harmonic
        this.J6 = 8.6e-5; // Sixth zonal harmonic
        
        // Magnetic field properties
        this.magneticField = {
            dipoleMoment: 4.6e18, // A⋅m² (Tesla⋅m³)
            tilt: THREE.MathUtils.degToRad(0.0), // Very small tilt compared to rotation axis
            offset: new THREE.Vector3(0, 0, 0), // Centered dipole
            strength: 0.21e-4, // Tesla at equator
            rotationPeriod: 38520 // Same as planet rotation
        };
        
        // Atmospheric properties
        this.atmosphere = {
            composition: {
                H2: 0.963, // Hydrogen
                He: 0.0325, // Helium  
                CH4: 4.5e-3, // Methane
                NH3: 1.25e-4, // Ammonia
                H2O: 1.1e-5, // Water vapor
                C2H6: 7e-6 // Ethane
            },
            scaleHeight: 60, // km
            surfacePressure: 101325, // Pa (1 bar at 1 bar level)
            temperatureProfile: {
                troposphere: 95, // K
                stratosphere: 160, // K
                thermosphere: 420 // K
            },
            windSpeeds: {
                equatorial: 500, // m/s (1800 km/h)
                temperate: 150, // m/s
                polar: 50 // m/s
            }
        };
        
        // Radiation and heat properties
        this.radiation = {
            solarConstant: 14.9, // W/m² (at Saturn's distance)
            internalHeat: 2.01, // Ratio of emitted to absorbed energy
            luminosity: 4.7e17, // W (total power radiated)
            effectiveTemperature: 95 // K
        };
        
        // Ring interaction properties
        this.ringSystem = null; // Will be set externally
        this.shepherdMoons = []; // Moons that interact with rings
        
        // Mesh components
        this.coreMesh = null;
        this.atmosphereMesh = null;
        this.cloudLayers = [];
        this.magnetosphereMesh = null;
        
        // Physics state
        this.atmosphericDynamics = {
            jetStreams: [],
            storms: [],
            hexagonalWave: null // North polar hexagon
        };
        
        // Initialize the planet
        this.initializePlanet();
    }
    
    /**
     * Initialize Saturn's visual representation and physics
     */
    initializePlanet() {
        this.createCoreMesh();
        this.createAtmosphere();
        this.createCloudLayers();
        this.initializeJetStreams();
        this.createMagnetosphere();
    }
    
    /**
     * Create Saturn's core mesh with oblate shape
     */
    createCoreMesh() {
        // Create oblate spheroid geometry
        const geometry = new THREE.SphereGeometry(1, 64, 64);
        
        // Modify vertices to create oblate shape
        const positions = geometry.attributes.position;
        for (let i = 0; i < positions.count; i++) {
            const y = positions.getY(i);
            const scaleFactor = this.polarRadius / this.radius;
            positions.setY(i, y * scaleFactor);
        }
        positions.needsUpdate = true;
        geometry.computeVertexNormals();
        
        // Create procedural surface texture
        const surfaceTexture = this.createSurfaceTexture();
        
        const material = new THREE.MeshPhongMaterial({
            map: surfaceTexture,
            emissive: new THREE.Color(0x332211),
            emissiveIntensity: 0.1,
            shininess: 10,
            transparent: false
        });
        
        this.coreMesh = new THREE.Mesh(geometry, material);
        this.coreMesh.castShadow = true;
        this.coreMesh.receiveShadow = true;
        this.coreMesh.rotation.z = this.axialTilt;
        
        this.mesh = this.coreMesh; // Primary mesh for base class
    }
    
    /**
     * Create procedural Saturn surface texture with bands
     */
    createSurfaceTexture() {
        const canvas = document.createElement('canvas');
        canvas.width = 2048;
        canvas.height = 1024;
        const ctx = canvas.getContext('2d');
        
        // Base color gradient
        const gradient = ctx.createLinearGradient(0, 0, 0, 1024);
        gradient.addColorStop(0, '#F5E6B8'); // North pole
        gradient.addColorStop(0.1, '#EDD9A3');
        gradient.addColorStop(0.2, '#E4C48E');
        gradient.addColorStop(0.3, '#D8B882');
        gradient.addColorStop(0.4, '#CEB077');
        gradient.addColorStop(0.5, '#C4A76D'); // Equator
        gradient.addColorStop(0.6, '#BB9E63');
        gradient.addColorStop(0.7, '#B39659');
        gradient.addColorStop(0.8, '#AA8D4F');
        gradient.addColorStop(0.9, '#A18545');
        gradient.addColorStop(1, '#987C3B'); // South pole
        
        ctx.fillStyle = gradient;
        ctx.fillRect(0, 0, 2048, 1024);
        
        // Add atmospheric bands
        this.addAtmosphericBands(ctx, 2048, 1024);
        
        // Add storm features
        this.addStormFeatures(ctx, 2048, 1024);
        
        // Add subtle noise for texture
        this.addNoiseTexture(ctx, 2048, 1024);
        
        return new THREE.CanvasTexture(canvas);
    }
    
    /**
     * Add atmospheric bands to texture
     */
    addAtmosphericBands(ctx, width, height) {
        const bandCount = 25;
        for (let i = 0; i < bandCount; i++) {
            const y = (i / bandCount) * height;
            const bandHeight = height / bandCount;
            
            // Vary band intensity
            const intensity = 0.05 + 0.03 * Math.sin(i * 0.7);
            const alpha = intensity * (0.8 + 0.4 * Math.random());
            
            // Alternate band colors
            const color = i % 2 === 0 ? `rgba(200, 180, 140, ${alpha})` : `rgba(180, 160, 120, ${alpha})`;
            
            ctx.fillStyle = color;
            
            // Create wavy band edges
            ctx.beginPath();
            ctx.moveTo(0, y);
            for (let x = 0; x <= width; x += 10) {
                const waveY = y + 3 * Math.sin(x * 0.01 + i * 0.5);
                ctx.lineTo(x, waveY);
            }
            ctx.lineTo(width, y + bandHeight);
            ctx.lineTo(0, y + bandHeight);
            ctx.closePath();
            ctx.fill();
        }
    }
    
    /**
     * Add storm features like the Great White Spot
     */
    addStormFeatures(ctx, width, height) {
        // Great White Spot (periodic storm)
        const stormX = width * 0.3;
        const stormY = height * 0.4;
        const stormWidth = width * 0.15;
        const stormHeight = height * 0.08;
        
        const gradient = ctx.createRadialGradient(
            stormX, stormY, 0,
            stormX, stormY, stormWidth
        );
        gradient.addColorStop(0, 'rgba(255, 245, 225, 0.6)');
        gradient.addColorStop(0.5, 'rgba(240, 220, 180, 0.3)');
        gradient.addColorStop(1, 'rgba(220, 200, 160, 0.1)');
        
        ctx.fillStyle = gradient;
        ctx.beginPath();
        ctx.ellipse(stormX, stormY, stormWidth, stormHeight, 0, 0, Math.PI * 2);
        ctx.fill();
        
        // Add smaller storms
        for (let i = 0; i < 20; i++) {
            const x = Math.random() * width;
            const y = Math.random() * height;
            const size = 10 + Math.random() * 30;
            const opacity = 0.1 + Math.random() * 0.2;
            
            ctx.fillStyle = `rgba(230, 210, 170, ${opacity})`;
            ctx.beginPath();
            ctx.ellipse(x, y, size * 2, size, Math.random() * Math.PI, 0, Math.PI * 2);
            ctx.fill();
        }
    }
    
    /**
     * Add noise texture for realism
     */
    addNoiseTexture(ctx, width, height) {
        const imageData = ctx.getImageData(0, 0, width, height);
        const data = imageData.data;
        
        for (let i = 0; i < data.length; i += 4) {
            const noise = (Math.random() - 0.5) * 10;
            data[i] = Math.max(0, Math.min(255, data[i] + noise));     // R
            data[i + 1] = Math.max(0, Math.min(255, data[i + 1] + noise)); // G
            data[i + 2] = Math.max(0, Math.min(255, data[i + 2] + noise)); // B
        }
        
        ctx.putImageData(imageData, 0, 0);
    }
    
    /**
     * Create atmospheric layers
     */
    createAtmosphere() {
        // Outer atmosphere
        const atmGeometry = new THREE.SphereGeometry(1.05, 64, 64);
        const atmMaterial = new THREE.MeshPhongMaterial({
            color: 0xFFE4B5,
            transparent: true,
            opacity: 0.15,
            emissive: 0xFFD700,
            emissiveIntensity: 0.02,
            side: THREE.BackSide
        });
        
        this.atmosphereMesh = new THREE.Mesh(atmGeometry, atmMaterial);
        this.atmosphereMesh.rotation.z = this.axialTilt;
    }
    
    /**
     * Create multiple cloud layers with different properties
     */
    createCloudLayers() {
        const layers = [
            { radius: 1.02, color: 0xF5E6B8, opacity: 0.3, speed: 1.0 },
            { radius: 1.03, color: 0xE4C48E, opacity: 0.25, speed: 0.8 },
            { radius: 1.04, color: 0xD8B882, opacity: 0.2, speed: 0.6 }
        ];
        
        layers.forEach((layer, index) => {
            const geometry = new THREE.SphereGeometry(layer.radius, 64, 64);
            const material = new THREE.MeshPhongMaterial({
                color: layer.color,
                transparent: true,
                opacity: layer.opacity,
                emissive: layer.color,
                emissiveIntensity: 0.01
            });
            
            const cloudMesh = new THREE.Mesh(geometry, material);
            cloudMesh.rotation.z = this.axialTilt;
            cloudMesh.userData = { rotationSpeed: layer.speed };
            
            this.cloudLayers.push(cloudMesh);
        });
    }
    
    /**
     * Initialize jet stream patterns
     */
    initializeJetStreams() {
        // Saturn has ~10-20 alternating jet streams
        const jetCount = 16;
        for (let i = 0; i < jetCount; i++) {
            const latitude = -80 + (i / (jetCount - 1)) * 160; // -80° to +80°
            const speed = this.calculateJetStreamSpeed(latitude);
            const direction = i % 2 === 0 ? 1 : -1; // Alternating directions
            
            this.atmosphericDynamics.jetStreams.push({
                latitude: THREE.MathUtils.degToRad(latitude),
                speed: speed * direction,
                width: THREE.MathUtils.degToRad(10), // 10 degree width
                strength: 1.0
            });
        }
        
        // Add equatorial superrotation
        this.atmosphericDynamics.jetStreams.push({
            latitude: 0,
            speed: 500, // 500 m/s eastward
            width: THREE.MathUtils.degToRad(15),
            strength: 1.5
        });
    }
    
    /**
     * Calculate jet stream speed based on latitude
     */
    calculateJetStreamSpeed(latitude) {
        const absLat = Math.abs(latitude);
        
        if (absLat < 10) return 500; // Equatorial jets
        if (absLat < 30) return 200; // Temperate jets
        if (absLat < 60) return 100; // Mid-latitude jets
        return 50; // Polar jets
    }
    
    /**
     * Create magnetosphere visualization
     */
    createMagnetosphere() {
        // Simplified magnetosphere shape
        const fieldLines = [];
        const lineCount = 24;
        
        for (let i = 0; i < lineCount; i++) {
            const phi = (i / lineCount) * Math.PI * 2;
            const points = [];
            
            // Create dipole field line
            for (let j = 0; j <= 50; j++) {
                const t = j / 50;
                const L = 3 + t * 8; // L-shell parameter
                const theta = Math.acos(Math.sqrt(t));
                
                const r = L * Math.sin(theta) * Math.sin(theta);
                const x = r * Math.cos(phi) * 3;
                const y = r * Math.sin(phi) * 3;
                const z = r * (Math.cos(theta) - Math.cos(Math.PI - theta)) * 1.5;
                
                points.push(new THREE.Vector3(x, z, y));
            }
            
            const geometry = new THREE.BufferGeometry().setFromPoints(points);
            const material = new THREE.LineBasicMaterial({
                color: 0x4466FF,
                transparent: true,
                opacity: 0.3
            });
            
            fieldLines.push(new THREE.Line(geometry, material));
        }
        
        this.magnetosphereMesh = new THREE.Group();
        fieldLines.forEach(line => this.magnetosphereMesh.add(line));
        this.magnetosphereMesh.visible = false; // Hidden by default
    }
    
    /**
     * Calculate gravitational acceleration including oblateness (J2, J4 effects)
     */
    calculateOblateness(position) {
        const r = position.length() * 1000; // Convert to meters
        const lat = Math.asin(position.y / (position.length())); // Latitude
        
        const sinLat = Math.sin(lat);
        const cosLat = Math.cos(lat);
        const P2 = 0.5 * (3 * sinLat * sinLat - 1); // Legendre polynomial P2
        const P4 = 0.125 * (35 * Math.pow(sinLat, 4) - 30 * sinLat * sinLat + 3); // P4
        
        const Re = this.radius * 1000; // Equatorial radius in meters
        const mu = this.gravParameter;
        
        // Radial component
        const ar = -mu / (r * r) * (1 + this.J2 * Math.pow(Re/r, 2) * P2 + this.J4 * Math.pow(Re/r, 4) * P4);
        
        // Latitudinal component  
        const dP2 = 3 * sinLat * cosLat;
        const dP4 = 2.5 * sinLat * cosLat * (7 * sinLat * sinLat - 3);
        const alat = -mu / (r * r) * sinLat * cosLat * (this.J2 * Math.pow(Re/r, 2) * dP2 + this.J4 * Math.pow(Re/r, 4) * dP4);
        
        // Convert back to Cartesian coordinates
        const unitR = position.clone().normalize();
        const unitLat = new THREE.Vector3(-unitR.x * unitR.y, 1 - unitR.y * unitR.y, -unitR.z * unitR.y).normalize();
        
        return unitR.multiplyScalar(ar / 1000).add(unitLat.multiplyScalar(alat / 1000)); // Convert to km/s²
    }
    
    /**
     * Update atmospheric dynamics
     */
    updateAtmosphere(deltaTime) {
        // Update cloud layer rotations
        this.cloudLayers.forEach(cloud => {
            cloud.rotation.y += deltaTime * cloud.userData.rotationSpeed * 0.05;
        });
        
        // Update jet streams (simplified)
        this.atmosphericDynamics.jetStreams.forEach(jet => {
            jet.strength *= (1 + 0.001 * (Math.random() - 0.5)); // Small fluctuations
        });
        
        // Rotate atmosphere mesh
        if (this.atmosphereMesh) {
            this.atmosphereMesh.rotation.y += deltaTime * 0.03; // Slightly slower than core
        }
    }
    
    /**
     * Update magnetosphere rotation
     */
    updateMagnetosphere(deltaTime) {
        if (this.magnetosphereMesh) {
            this.magnetosphereMesh.rotation.y += deltaTime * (2 * Math.PI / this.magneticField.rotationPeriod);
        }
    }
    
    /**
     * Get atmospheric density at altitude
     */
    getAtmosphericDensity(altitude) {
        // Exponential atmosphere model
        const scaleHeight = this.atmosphere.scaleHeight * 1000; // Convert to meters
        const surfaceDensity = 0.19; // kg/m³ at 1 bar level
        return surfaceDensity * Math.exp(-altitude / scaleHeight);
    }
    
    /**
     * Get temperature at altitude
     */
    getTemperatureAtAltitude(altitude) {
        if (altitude < 50000) {
            // Troposphere
            return this.atmosphere.temperatureProfile.troposphere - altitude * 0.002;
        } else if (altitude < 200000) {
            // Stratosphere  
            return this.atmosphere.temperatureProfile.stratosphere;
        } else {
            // Thermosphere
            return this.atmosphere.temperatureProfile.thermosphere;
        }
    }
    
    /**
     * Calculate magnetic field at position
     */
    getMagneticField(position) {
        const r = position.length() * 1000; // Convert to meters
        const Re = this.radius * 1000;
        
        if (r < Re) return new THREE.Vector3(0, 0, 0); // Inside planet
        
        // Dipole field approximation
        const m = this.magneticField.dipoleMoment;
        const mu0 = 4 * Math.PI * 1e-7; // Permeability of free space
        
        const unitR = position.clone().normalize();
        const dipoleAxis = new THREE.Vector3(0, 1, 0); // Aligned with rotation axis
        
        const cosTheta = unitR.dot(dipoleAxis);
        const Br = (mu0 * m / (4 * Math.PI)) * (2 * cosTheta) / Math.pow(r, 3);
        const Btheta = (mu0 * m / (4 * Math.PI)) * Math.sin(Math.acos(cosTheta)) / Math.pow(r, 3);
        
        // Convert to Cartesian
        const sinTheta = Math.sqrt(1 - cosTheta * cosTheta);
        const Bx = Br * sinTheta * unitR.x + Btheta * cosTheta * unitR.x;
        const By = Br * cosTheta;
        const Bz = Br * sinTheta * unitR.z + Btheta * cosTheta * unitR.z;
        
        return new THREE.Vector3(Bx, By, Bz);
    }
    
    /**
     * Override update method to include Saturn-specific physics
     */
    update(deltaTime) {
        // Call parent update
        super.updateRotation(deltaTime);
        
        // Update Saturn-specific systems
        this.updateAtmosphere(deltaTime);
        this.updateMagnetosphere(deltaTime);
        
        // Update internal heat generation
        this.updateInternalHeat(deltaTime);
    }
    
    /**
     * Update internal heat generation
     */
    updateInternalHeat(deltaTime) {
        // Saturn radiates 2.01 times more energy than it receives from the Sun
        // This is due to gravitational contraction and helium rain
        const solarInput = this.radiation.solarConstant * Math.PI * Math.pow(this.radius * 1000, 2);
        const internalInput = solarInput * (this.radiation.internalHeat - 1);
        
        // Simple temperature regulation (simplified)
        const targetTemp = this.radiation.effectiveTemperature;
        this.temperature += (targetTemp - this.temperature) * deltaTime * 0.001;
    }
    
    /**
     * Add Saturn and its components to scene
     */
    addToScene(scene) {
        if (this.coreMesh) scene.add(this.coreMesh);
        if (this.atmosphereMesh) scene.add(this.atmosphereMesh);
        this.cloudLayers.forEach(cloud => scene.add(cloud));
        if (this.magnetosphereMesh) scene.add(this.magnetosphereMesh);
    }
    
    /**
     * Remove Saturn from scene
     */
    removeFromScene(scene) {
        if (this.coreMesh) scene.remove(this.coreMesh);
        if (this.atmosphereMesh) scene.remove(this.atmosphereMesh);
        this.cloudLayers.forEach(cloud => scene.remove(cloud));
        if (this.magnetosphereMesh) scene.remove(this.magnetosphereMesh);
    }
    
    /**
     * Toggle magnetosphere visibility
     */
    toggleMagnetosphere(visible) {
        if (this.magnetosphereMesh) {
            this.magnetosphereMesh.visible = visible;
        }
    }
    
    /**
     * Toggle atmosphere visibility
     */
    toggleAtmosphere(visible) {
        if (this.atmosphereMesh) {
            this.atmosphereMesh.visible = visible;
        }
        this.cloudLayers.forEach(cloud => {
            cloud.visible = visible;
        });
    }
    
    /**
     * Get detailed planet information
     */
    getInfo() {
        return {
            ...super.exportState(),
            oblateness: (this.radius - this.polarRadius) / this.radius,
            gravitationalHarmonics: { J2: this.J2, J4: this.J4, J6: this.J6 },
            magneticField: this.magneticField,
            atmosphere: this.atmosphere,
            radiation: this.radiation,
            jetStreamCount: this.atmosphericDynamics.jetStreams.length,
            internalHeat: this.radiation.internalHeat
        };
    }
    
    /**
     * Dispose of all resources
     */
    dispose() {
        super.dispose();
        
        if (this.atmosphereMesh) {
            this.atmosphereMesh.geometry.dispose();
            this.atmosphereMesh.material.dispose();
        }
        
        this.cloudLayers.forEach(cloud => {
            cloud.geometry.dispose();
            cloud.material.dispose();
        });
        
        if (this.magnetosphereMesh) {
            this.magnetosphereMesh.children.forEach(child => {
                child.geometry.dispose();
                child.material.dispose();
            });
        }
    }
}
