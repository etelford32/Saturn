/**
 * Engine.js - Main Saturn Simulation Engine
 * Coordinates all subsystems, manages time stepping, and ensures deterministic execution
 */

import * as THREE from 'three';

export class Engine {
    constructor(config = {}) {
        // Core engine configuration
        this.config = {
            targetFPS: config.targetFPS || 60,
            maxDeltaTime: config.maxDeltaTime || 1/30, // Cap at 30fps minimum
            timescale: config.timescale || 1.0,
            enablePhysics: config.enablePhysics !== false,
            enableDeterministic: config.enableDeterministic !== false,
            qualityLevel: config.qualityLevel || 'medium',
            enableDebug: config.enableDebug || false,
            ...config
        };
        
        // Time management
        this.time = {
            current: 0,           // Current simulation time (seconds)
            delta: 0,             // Time since last frame (seconds)
            scale: this.config.timescale,
            accumulator: 0,       // For fixed timestep physics
            alpha: 0,             // Interpolation factor for rendering
            frameCount: 0,
            realTime: 0,          // Real world time elapsed
            startTime: performance.now()
        };
        
        // Multi-rate time stepping for different subsystems
        this.timeSteps = {
            orbital: 0.1,         // Orbital mechanics (100ms)
            rings: 0.01,          // Ring particle physics (10ms)  
            atmosphere: 1.0,      // Atmospheric dynamics (1s)
            magnetosphere: 0.5,   // Magnetic field updates (500ms)
            thermal: 10.0,        // Thermal evolution (10s)
            render: 1/60          // Rendering rate (60fps)
        };
        
        // Subsystem management
        this.subsystems = new Map();
        this.subsystemOrder = []; // Execution order
        this.subsystemStates = new Map(); // Track update timing
        
        // Registered objects and their update schedules
        this.celestialBodies = [];
        this.ringSystem = null;
        this.atmosphereSystem = null;
        this.magnetosphereSystem = null;
        
        // Physics integration
        this.physics = {
            integrator: config.integrator || 'leapfrog', // leapfrog, rk4, verlet
            energyConservation: true,
            energyDrift: 0,
            energyDriftThreshold: 1e-6,
            lastTotalEnergy: 0,
            collisionDetection: config.enableCollisions !== false,
            selfGravity: config.enableSelfGravity !== false
        };
        
        // Performance monitoring
        this.performance = {
            frameTime: 0,
            physicsTime: 0,
            renderTime: 0,
            subsystemTimes: new Map(),
            averageFrameTime: 0,
            frameTimeHistory: [],
            maxHistoryLength: 100,
            bottlenecks: [],
            memoryUsage: 0
        };
        
        // State management
        this.state = {
            running: false,
            paused: false,
            initialized: false,
            error: null,
            saveStateInterval: 60, // seconds
            lastSaveTime: 0,
            autoSave: config.autoSave !== false
        };
        
        // Deterministic execution
        this.deterministic = {
            seed: config.seed || Date.now(),
            rng: null, // Will be initialized
            frameSeeds: new Map(), // Per-frame seeds for each subsystem
            reproducible: this.config.enableDeterministic
        };
        
        // Event system
        this.events = {
            listeners: new Map(),
            queue: [],
            maxQueueSize: 1000
        };
        
        // Debug and diagnostics
        this.debug = {
            enabled: this.config.enableDebug,
            logLevel: config.logLevel || 'info',
            frameGraph: [],
            energyGraph: [],
            performanceGraph: [],
            memoryGraph: []
        };
        
        // Quality level settings
        this.qualitySettings = this.getQualitySettings(this.config.qualityLevel);
        
        this.initialize();
    }
    
    /**
     * Initialize the simulation engine
     */
    initialize() {
        this.initializeDeterministicRNG();
        this.initializeSubsystems();
        this.initializePerformanceMonitoring();
        this.setupEventHandlers();
        
        this.state.initialized = true;
        this.emit('engine:initialized');
        
        this.log('info', 'Saturn Simulation Engine initialized', {
            quality: this.config.qualityLevel,
            deterministic: this.deterministic.reproducible,
            seed: this.deterministic.seed
        });
    }
    
    /**
     * Initialize deterministic random number generator
     */
    initializeDeterministicRNG() {
        if (this.deterministic.reproducible) {
            // Simple LCG (Linear Congruential Generator) for deterministic randomness
            this.deterministic.rng = {
                seed: this.deterministic.seed,
                next: function() {
                    this.seed = (this.seed * 1664525 + 1013904223) % 4294967296;
                    return this.seed / 4294967296;
                }
            };
        }
    }
    
    /**
     * Get quality settings based on level
     */
    getQualitySettings(level) {
        const settings = {
            low: {
                ringParticles: 5000,
                atmosphereResolution: 64,
                physicsSubsteps: 1,
                enableSelfGravity: false,
                enableAtmosphere: false,
                enableMagnetosphere: false,
                collisionDetection: 'basic'
            },
            medium: {
                ringParticles: 25000,
                atmosphereResolution: 128,
                physicsSubsteps: 2,
                enableSelfGravity: true,
                enableAtmosphere: true,
                enableMagnetosphere: false,
                collisionDetection: 'spatial_hash'
            },
            high: {
                ringParticles: 75000,
                atmosphereResolution: 256,
                physicsSubsteps: 4,
                enableSelfGravity: true,
                enableAtmosphere: true,
                enableMagnetosphere: true,
                collisionDetection: 'spatial_hash'
            },
            ultra: {
                ringParticles: 200000,
                atmosphereResolution: 512,
                physicsSubsteps: 8,
                enableSelfGravity: true,
                enableAtmosphere: true,
                enableMagnetosphere: true,
                collisionDetection: 'gpu_accelerated'
            }
        };
        
        return settings[level] || settings.medium;
    }
    
    /**
     * Initialize subsystem management
     */
    initializeSubsystems() {
        // Define subsystem execution order
        this.subsystemOrder = [
            'orbital',        // Gravitational dynamics first
            'rings',          // Ring particle physics
            'atmosphere',     // Atmospheric dynamics
            'magnetosphere',  // Magnetic field interactions
            'thermal'         // Thermal evolution last
        ];
        
        // Initialize subsystem states
        this.subsystemOrder.forEach(name => {
            this.subsystemStates.set(name, {
                lastUpdate: 0,
                accumulator: 0,
                updateCount: 0,
                averageTime: 0,
                enabled: true
            });
        });
    }
    
    /**
     * Initialize performance monitoring
     */
    initializePerformanceMonitoring() {
        // Initialize performance timers for each subsystem
        this.subsystemOrder.forEach(name => {
            this.performance.subsystemTimes.set(name, {
                current: 0,
                average: 0,
                maximum: 0,
                history: []
            });
        });
        
        // Start memory monitoring if available
        if (performance.memory) {
            setInterval(() => {
                this.performance.memoryUsage = performance.memory.usedJSHeapSize;
            }, 1000);
        }
    }
    
    /**
     * Setup event handlers
     */
    setupEventHandlers() {
        // Handle window visibility changes
        if (typeof document !== 'undefined') {
            document.addEventListener('visibilitychange', () => {
                if (document.hidden) {
                    this.pause();
                } else {
                    this.resume();
                }
            });
        }
        
        // Handle errors
        window.addEventListener('error', (event) => {
            this.handleError(event.error);
        });
    }
    
    /**
     * Register a subsystem
     */
    registerSubsystem(name, subsystem, config = {}) {
        this.subsystems.set(name, {
            instance: subsystem,
            enabled: config.enabled !== false,
            updateRate: config.updateRate || this.timeSteps[name] || 1/60,
            priority: config.priority || 0,
            dependencies: config.dependencies || []
        });
        
        this.log('info', `Registered subsystem: ${name}`);
        this.emit('subsystem:registered', { name, subsystem });
    }
    
    /**
     * Register a celestial body
     */
    registerCelestialBody(body) {
        this.celestialBodies.push(body);
        this.log('debug', `Registered celestial body: ${body.name}`);
        this.emit('body:registered', { body });
    }
    
    /**
     * Register ring system
     */
    registerRingSystem(ringSystem) {
        this.ringSystem = ringSystem;
        this.registerSubsystem('rings', ringSystem, {
            updateRate: this.timeSteps.rings,
            priority: 2
        });
    }
    
    /**
     * Register atmosphere system
     */
    registerAtmosphereSystem(atmosphere) {
        this.atmosphereSystem = atmosphere;
        this.registerSubsystem('atmosphere', atmosphere, {
            updateRate: this.timeSteps.atmosphere,
            priority: 3
        });
    }
    
    /**
     * Start the simulation
     */
    start() {
        if (!this.state.initialized) {
            throw new Error('Engine not initialized');
        }
        
        this.state.running = true;
        this.state.paused = false;
        this.time.startTime = performance.now();
        
        this.log('info', 'Simulation started');
        this.emit('engine:started');
        
        this.animate();
    }
    
    /**
     * Stop the simulation
     */
    stop() {
        this.state.running = false;
        this.state.paused = false;
        
        this.log('info', 'Simulation stopped');
        this.emit('engine:stopped');
    }
    
    /**
     * Pause the simulation
     */
    pause() {
        this.state.paused = true;
        this.log('info', 'Simulation paused');
        this.emit('engine:paused');
    }
    
    /**
     * Resume the simulation
     */
    resume() {
        this.state.paused = false;
        this.log('info', 'Simulation resumed');
        this.emit('engine:resumed');
    }
    
    /**
     * Reset the simulation
     */
    reset() {
        this.stop();
        
        // Reset time
        this.time.current = 0;
        this.time.frameCount = 0;
        this.time.accumulator = 0;
        
        // Reset subsystem states
        this.subsystemStates.forEach(state => {
            state.lastUpdate = 0;
            state.accumulator = 0;
            state.updateCount = 0;
        });
        
        // Reset performance
        this.performance.frameTimeHistory = [];
        this.performance.subsystemTimes.forEach(timing => {
            timing.history = [];
            timing.average = 0;
            timing.maximum = 0;
        });
        
        // Reset celestial bodies
        this.celestialBodies.forEach(body => {
            if (body.reset) body.reset();
        });
        
        this.log('info', 'Simulation reset');
        this.emit('engine:reset');
    }
    
    /**
     * Main animation loop
     */
    animate() {
        if (!this.state.running) return;
        
        requestAnimationFrame(() => this.animate());
        
        const frameStart = performance.now();
        
        try {
            // Calculate frame timing
            this.updateTiming(frameStart);
            
            if (!this.state.paused && this.time.delta > 0) {
                // Update simulation
                this.updateSimulation();
                
                // Monitor energy conservation
                if (this.physics.energyConservation) {
                    this.monitorEnergyConservation();
                }
                
                // Auto-save state periodically
                if (this.state.autoSave && 
                    this.time.current - this.state.lastSaveTime > this.state.saveStateInterval) {
                    this.saveState();
                }
            }
            
            // Update performance metrics
            this.updatePerformanceMetrics(frameStart);
            
            // Process event queue
            this.processEventQueue();
            
        } catch (error) {
            this.handleError(error);
        }
    }
    
    /**
     * Update timing information
     */
    updateTiming(currentTime) {
        const realDelta = (currentTime - this.time.startTime) / 1000 - this.time.realTime;
        this.time.realTime = (currentTime - this.time.startTime) / 1000;
        
        // Cap delta time to prevent spiral of death
        this.time.delta = Math.min(realDelta * this.time.scale, this.config.maxDeltaTime);
        this.time.current += this.time.delta;
        this.time.frameCount++;
        
        // Update accumulators for fixed timestep systems
        this.subsystemStates.forEach((state, name) => {
            const subsystem = this.subsystems.get(name);
            if (subsystem) {
                state.accumulator += this.time.delta;
            }
        });
    }
    
    /**
     * Update the entire simulation
     */
    updateSimulation() {
        // Generate deterministic seeds for this frame
        if (this.deterministic.reproducible) {
            this.generateFrameSeeds();
        }
        
        // Update subsystems in order
        this.subsystemOrder.forEach(name => {
            this.updateSubsystem(name);
        });
        
        // Update celestial bodies (orbital mechanics)
        this.updateCelestialBodies();
        
        // Emit frame update event
        this.emit('engine:frame', {
            time: this.time.current,
            delta: this.time.delta,
            frame: this.time.frameCount
        });
    }
    
    /**
     * Update a specific subsystem with appropriate timestep
     */
    updateSubsystem(name) {
        const subsystem = this.subsystems.get(name);
        const state = this.subsystemStates.get(name);
        
        if (!subsystem || !subsystem.enabled || !state.enabled) return;
        
        const startTime = performance.now();
        
        try {
            const targetDt = subsystem.updateRate;
            
            // Fixed timestep with accumulator
            while (state.accumulator >= targetDt) {
                // Set deterministic seed for this subsystem
                if (this.deterministic.reproducible) {
                    const seed = this.deterministic.frameSeeds.get(name);
                    if (seed !== undefined) {
                        Math.seedrandom = () => seed;
                    }
                }
                
                // Update with fixed timestep
                if (subsystem.instance && subsystem.instance.update) {
                    subsystem.instance.update(targetDt, this.time.scale);
                }
                
                state.accumulator -= targetDt;
                state.lastUpdate = this.time.current;
                state.updateCount++;
                
                // Prevent infinite loops
                if (state.accumulator > targetDt * 10) {
                    state.accumulator = 0;
                    break;
                }
            }
            
            // Calculate interpolation factor for rendering
            this.time.alpha = state.accumulator / targetDt;
            
        } catch (error) {
            this.log('error', `Subsystem ${name} update failed:`, error);
            state.enabled = false; // Disable failing subsystem
        }
        
        // Record timing
        const updateTime = performance.now() - startTime;
        this.recordSubsystemTiming(name, updateTime);
    }
    
    /**
     * Update celestial bodies (orbital mechanics)
     */
    updateCelestialBodies() {
        const dt = this.timeSteps.orbital;
        
        // Update orbital mechanics for all bodies
        this.celestialBodies.forEach(body => {
            if (body.update) {
                body.update(dt * this.time.scale);
            }
        });
        
        // Calculate gravitational interactions
        if (this.physics.integrator && this.celestialBodies.length > 1) {
            this.updateGravitationalForces(dt);
        }
    }
    
    /**
     * Update gravitational forces between bodies
     */
    updateGravitationalForces(dt) {
        // Clear existing forces
        this.celestialBodies.forEach(body => {
            if (body.clearForces) body.clearForces();
        });
        
        // N-body gravitational interactions
        for (let i = 0; i < this.celestialBodies.length; i++) {
            for (let j = i + 1; j < this.celestialBodies.length; j++) {
                const bodyA = this.celestialBodies[i];
                const bodyB = this.celestialBodies[j];
                
                if (bodyA.calculateGravitationalAcceleration && bodyB.addForce) {
                    const forceA = bodyA.calculateGravitationalAcceleration(bodyB);
                    const forceB = bodyB.calculateGravitationalAcceleration(bodyA);
                    
                    bodyA.addForce(forceA);
                    bodyB.addForce(forceB);
                }
            }
        }
        
        // Apply physics integration
        this.celestialBodies.forEach(body => {
            if (body.updatePhysics) {
                body.updatePhysics(dt, this.physics.integrator);
            }
        });
    }
    
    /**
     * Generate deterministic seeds for each subsystem
     */
    generateFrameSeeds() {
        this.deterministic.frameSeeds.clear();
        
        let baseSeed = this.deterministic.seed + this.time.frameCount;
        
        this.subsystemOrder.forEach(name => {
            baseSeed = (baseSeed * 1664525 + 1013904223) % 4294967296;
            this.deterministic.frameSeeds.set(name, baseSeed / 4294967296);
        });
    }
    
    /**
     * Monitor energy conservation
     */
    monitorEnergyConservation() {
        let totalEnergy = 0;
        
        // Calculate total system energy
        this.celestialBodies.forEach(body => {
            if (body.velocity && body.mass) {
                // Kinetic energy
                const kineticEnergy = 0.5 * body.mass * body.velocity.lengthSq();
                totalEnergy += kineticEnergy;
                
                // Potential energy (simplified)
                if (body.position) {
                    const r = body.position.length() * 1000; // Convert to meters
                    if (r > 0 && body.parent && body.parent.gravParameter) {
                        const potentialEnergy = -body.parent.gravParameter * body.mass / r;
                        totalEnergy += potentialEnergy;
                    }
                }
            }
        });
        
        // Check energy drift
        if (this.physics.lastTotalEnergy !== 0) {
            this.physics.energyDrift = Math.abs(
                (totalEnergy - this.physics.lastTotalEnergy) / this.physics.lastTotalEnergy
            );
            
            if (this.physics.energyDrift > this.physics.energyDriftThreshold) {
                this.log('warn', `Energy conservation violated: drift = ${this.physics.energyDrift}`);
                this.emit('physics:energy_drift', { drift: this.physics.energyDrift });
            }
        }
        
        this.physics.lastTotalEnergy = totalEnergy;
        
        // Record for debugging
        if (this.debug.enabled) {
            this.debug.energyGraph.push({
                time: this.time.current,
                energy: totalEnergy,
                drift: this.physics.energyDrift
            });
            
            // Limit history size
            if (this.debug.energyGraph.length > 1000) {
                this.debug.energyGraph.shift();
            }
        }
    }
    
    /**
     * Record subsystem timing
     */
    recordSubsystemTiming(name, time) {
        const timing = this.performance.subsystemTimes.get(name);
        if (timing) {
            timing.current = time;
            timing.history.push(time);
            
            // Limit history
            if (timing.history.length > 100) {
                timing.history.shift();
            }
            
            // Update statistics
            timing.average = timing.history.reduce((a, b) => a + b, 0) / timing.history.length;
            timing.maximum = Math.max(timing.maximum, time);
        }
    }
    
    /**
     * Update performance metrics
     */
    updatePerformanceMetrics(frameStart) {
        this.performance.frameTime = performance.now() - frameStart;
        this.performance.frameTimeHistory.push(this.performance.frameTime);
        
        // Limit history
        if (this.performance.frameTimeHistory.length > this.performance.maxHistoryLength) {
            this.performance.frameTimeHistory.shift();
        }
        
        // Calculate average
        this.performance.averageFrameTime = 
            this.performance.frameTimeHistory.reduce((a, b) => a + b, 0) / 
            this.performance.frameTimeHistory.length;
        
        // Detect bottlenecks
        this.detectBottlenecks();
        
        // Record for debug graphs
        if (this.debug.enabled) {
            this.debug.performanceGraph.push({
                time: this.time.current,
                frameTime: this.performance.frameTime,
                physicsTime: this.performance.physicsTime
            });
            
            if (this.debug.performanceGraph.length > 1000) {
                this.debug.performanceGraph.shift();
            }
        }
    }
    
    /**
     * Detect performance bottlenecks
     */
    detectBottlenecks() {
        const threshold = 16.67; // 60fps threshold in ms
        
        this.performance.bottlenecks = [];
        
        // Check overall frame time
        if (this.performance.frameTime > threshold) {
            this.performance.bottlenecks.push({
                type: 'frame',
                severity: this.performance.frameTime / threshold
            });
        }
        
        // Check subsystem times
        this.performance.subsystemTimes.forEach((timing, name) => {
            if (timing.current > threshold * 0.25) { // 25% of frame budget
                this.performance.bottlenecks.push({
                    type: 'subsystem',
                    name: name,
                    severity: timing.current / (threshold * 0.25)
                });
            }
        });
    }
    
    /**
     * Save simulation state
     */
    saveState() {
        try {
            const state = {
                time: this.time,
                celestialBodies: this.celestialBodies.map(body => 
                    body.exportState ? body.exportState() : null
                ).filter(s => s !== null),
                ringSystem: this.ringSystem ? this.ringSystem.exportState() : null,
                atmosphere: this.atmosphereSystem ? this.atmosphereSystem.exportState() : null,
                metadata: {
                    version: '1.0.0',
                    timestamp: Date.now(),
                    seed: this.deterministic.seed
                }
            };
            
            // Save to localStorage or emit event for external handling
            if (typeof localStorage !== 'undefined') {
                localStorage.setItem('saturn_simulation_state', JSON.stringify(state));
            }
            
            this.state.lastSaveTime = this.time.current;
            this.emit('state:saved', { state });
            
        } catch (error) {
            this.log('error', 'Failed to save state:', error);
        }
    }
    
    /**
     * Load simulation state
     */
    loadState(state) {
        try {
            if (typeof state === 'string') {
                state = JSON.parse(state);
            }
            
            // Restore time
            if (state.time) {
                this.time = { ...this.time, ...state.time };
            }
            
            // Restore celestial bodies
            if (state.celestialBodies) {
                state.celestialBodies.forEach((bodyState, index) => {
                    if (this.celestialBodies[index] && this.celestialBodies[index].importState) {
                        this.celestialBodies[index].importState(bodyState);
                    }
                });
            }
            
            // Restore subsystems
            if (state.ringSystem && this.ringSystem && this.ringSystem.importState) {
                this.ringSystem.importState(state.ringSystem);
            }
            
            if (state.atmosphere && this.atmosphereSystem && this.atmosphereSystem.importState) {
                this.atmosphereSystem.importState(state.atmosphere);
            }
            
            this.emit('state:loaded', { state });
            this.log('info', 'State loaded successfully');
            
        } catch (error) {
            this.log('error', 'Failed to load state:', error);
            throw error;
        }
    }
    
    /**
     * Event system - emit event
     */
    emit(eventName, data = {}) {
        const listeners = this.events.listeners.get(eventName);
        if (listeners) {
            listeners.forEach(callback => {
                try {
                    callback(data);
                } catch (error) {
                    this.log('error', `Event listener error for ${eventName}:`, error);
                }
            });
        }
    }
    
    /**
     * Event system - add listener
     */
    on(eventName, callback) {
        if (!this.events.listeners.has(eventName)) {
            this.events.listeners.set(eventName, []);
        }
        this.events.listeners.get(eventName).push(callback);
    }
    
    /**
     * Event system - remove listener
     */
    off(eventName, callback) {
        const listeners = this.events.listeners.get(eventName);
        if (listeners) {
            const index = listeners.indexOf(callback);
            if (index > -1) {
                listeners.splice(index, 1);
            }
        }
    }
    
    /**
     * Process event queue
     */
    processEventQueue() {
        while (this.events.queue.length > 0) {
            const event = this.events.queue.shift();
            this.emit(event.name, event.data);
        }
    }
    
    /**
     * Handle errors
     */
    handleError(error) {
        this.state.error = error;
        this.log('error', 'Engine error:', error);
        this.emit('engine:error', { error });
        
        // Attempt graceful recovery
        if (this.state.running) {
            this.pause();
        }
    }
    
    /**
     * Logging system
     */
    log(level, message, data = null) {
        if (!this.debug.enabled && level === 'debug') return;
        
        const timestamp = new Date().toISOString();
        const logEntry = {
            timestamp,
            level,
            message,
            data,
            frame: this.time.frameCount,
            time: this.time.current
        };
        
        // Console output
        console[level] || console.log(`[${timestamp}] ${level.toUpperCase()}: ${message}`, data || '');
        
        // Emit log event for external handling
        this.emit('engine:log', logEntry);
    }
    
    /**
     * Get engine statistics
     */
    getStatistics() {
        return {
            time: this.time,
            performance: {
                ...this.performance,
                fps: 1000 / this.performance.averageFrameTime,
                targetFPS: this.config.targetFPS
            },
            physics: this.physics,
            subsystems: Object.fromEntries(
                Array.from(this.subsystemStates.entries()).map(([name, state]) => [
                    name, 
                    {
                        ...state,
                        timing: this.performance.subsystemTimes.get(name)
                    }
                ])
            ),
            memory: this.performance.memoryUsage,
            quality: this.config.qualityLevel,
            deterministic: this.deterministic.reproducible
        };
    }
    
    /**
     * Update quality level
     */
    setQualityLevel(level) {
        this.config.qualityLevel = level;
        this.qualitySettings = this.getQualitySettings(level);
        
        // Notify subsystems of quality change
        this.emit('engine:quality_changed', { 
            level, 
            settings: this.qualitySettings 
        });
        
        this.log('info', `Quality level changed to: ${level}`);
    }
    
    /**
     * Update time scale
     */
    setTimeScale(scale) {
        this.time.scale = Math.max(0, scale);
        this.emit('engine:timescale_changed', { scale: this.time.scale });
        this.log('info', `Time scale changed to: ${scale}x`);
    }
    
    /**
     * Get debug information
     */
    getDebugInfo() {
        return {
            ...this.debug,
            state: this.state,
            config: this.config,
            subsystems: Array.from(this.subsystems.keys()),
            celestialBodies: this.celestialBodies.map(body => body.name),
            bottlenecks: this.performance.bottlenecks
        };
    }
    
    /**
     * Dispose engine and clean up resources
     */
    dispose() {
        this.stop();
        
        // Clear all listeners
        this.events.listeners.clear();
        this.events.queue = [];
        
        // Dispose subsystems
        this.subsystems.forEach(subsystem => {
            if (subsystem.instance && subsystem.instance.dispose) {
                subsystem.instance.dispose();
            }
        });
        
        // Dispose celestial bodies
        this.celestialBodies.forEach(body => {
            if (body.dispose) {
                body.dispose();
            }
        });
        
        this.log('info', 'Engine disposed');
        this.emit('engine:disposed');
    }
}
