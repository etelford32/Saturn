/**
 * Saturn Simulation - Main Entry Point
 * Initializes Three.js scene and connects UI components
 */

import * as THREE from 'three';
import Stats from 'stats.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';

// Import simulation modules (these will be created next)
// import { Engine } from './core/Engine.js';
// import { Saturn } from './celestial/Saturn.js';
// import { RingSystem } from './celestial/PlanetRings.js';
// import { MoonSystem } from './celestial/MoonSystem.js';
// import { CameraController } from './camera/CameraController.js';
// import { UIManager } from './ui/UIManager.js';

class SaturnSimulation {
    constructor() {
        this.container = document.getElementById('canvas');
        this.loadingScreen = document.getElementById('loading-screen');
        this.loadingProgress = document.querySelector('.loading-progress-bar');
        this.loadingStatus = document.querySelector('.loading-status');
        
        // Simulation state
        this.isRunning = false;
        this.isPaused = false;
        this.timeScale = 1.0;
        this.elapsedTime = 0;
        this.deltaTime = 0;
        this.lastTime = 0;
        
        // Three.js components
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.composer = null;
        this.controls = null;
        this.stats = null;
        
        // Simulation objects
        this.saturn = null;
        this.ringSystem = null;
        this.moons = [];
        this.starfield = null;
        
        // Settings
        this.settings = {
            showRings: true,
            showMoons: true,
            showAtmosphere: true,
            showMagnetosphere: false,
            showOrbits: false,
            showLabels: true,
            showGrid: false,
            showStats: false,
            quality: 'medium'
        };
        
        this.init();
    }
    
    async init() {
        try {
            this.updateLoadingStatus('Initializing renderer...');
            await this.initRenderer();
            this.updateLoadingProgress(20);
            
            this.updateLoadingStatus('Creating scene...');
            await this.initScene();
            this.updateLoadingProgress(40);
            
            this.updateLoadingStatus('Loading Saturn...');
            await this.createSaturn();
            this.updateLoadingProgress(60);
            
            this.updateLoadingStatus('Creating ring system...');
            await this.createRings();
            this.updateLoadingProgress(70);
            
            this.updateLoadingStatus('Placing moons...');
            await this.createMoons();
            this.updateLoadingProgress(80);
            
            this.updateLoadingStatus('Setting up controls...');
            this.initControls();
            this.updateLoadingProgress(90);
            
            this.updateLoadingStatus('Initializing UI...');
            this.initUI();
            this.updateLoadingProgress(100);
            
            // Hide loading screen
            setTimeout(() => {
                this.hideLoadingScreen();
                this.start();
            }, 500);
            
        } catch (error) {
            console.error('Failed to initialize simulation:', error);
            this.updateLoadingStatus('Failed to initialize. Please refresh.');
        }
    }
    
    async initRenderer() {
        // Create renderer
        this.renderer = new THREE.WebGLRenderer({
            canvas: this.container,
            antialias: true,
            powerPreference: 'high-performance'
        });
        
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
        this.renderer.toneMappingExposure = 1.2;
        this.renderer.shadowMap.enabled = true;
        this.renderer.shadowMap.type = THREE.PCFSoftShadowMap;
        
        // Handle window resize
        window.addEventListener('resize', () => this.onWindowResize());
    }
    
    async initScene() {
        // Create scene
        this.scene = new THREE.Scene();
        this.scene.fog = new THREE.Fog(0x000011, 50, 500);
        
        // Create camera
        const aspect = window.innerWidth / window.innerHeight;
        this.camera = new THREE.PerspectiveCamera(60, aspect, 0.1, 1000);
        this.camera.position.set(30, 15, 30);
        this.camera.lookAt(0, 0, 0);
        
        // Add lights
        this.setupLighting();
        
        // Create starfield
        this.createStarfield();
        
        // Setup post-processing
        this.setupPostProcessing();
        
        // Setup performance stats
        if (typeof Stats !== 'undefined') {
            this.stats = new Stats();
            this.stats.showPanel(0);
            this.stats.dom.style.position = 'absolute';
            this.stats.dom.style.top = '10px';
            this.stats.dom.style.left = '10px';
            this.stats.dom.style.display = 'none';
            document.body.appendChild(this.stats.dom);
        }
    }
    
    setupLighting() {
        // Ambient light
        const ambientLight = new THREE.AmbientLight(0x111133, 0.4);
        this.scene.add(ambientLight);
        
        // Main sun light
        const sunLight = new THREE.DirectionalLight(0xffffff, 2);
        sunLight.position.set(100, 50, 50);
        sunLight.castShadow = true;
        sunLight.shadow.mapSize.width = 2048;
        sunLight.shadow.mapSize.height = 2048;
        sunLight.shadow.camera.near = 0.5;
        sunLight.shadow.camera.far = 500;
        sunLight.shadow.camera.left = -50;
        sunLight.shadow.camera.right = 50;
        sunLight.shadow.camera.top = 50;
        sunLight.shadow.camera.bottom = -50;
        this.scene.add(sunLight);
        
        // Saturn glow light
        const saturnGlow = new THREE.PointLight(0xffdd88, 0.5, 100);
        saturnGlow.position.set(0, 0, 0);
        this.scene.add(saturnGlow);
        
        // Rim light for better visibility
        const rimLight = new THREE.DirectionalLight(0x4466ff, 0.3);
        rimLight.position.set(-50, 0, -50);
        this.scene.add(rimLight);
    }
    
    createStarfield() {
        const starsGeometry = new THREE.BufferGeometry();
        const starCount = 15000;
        const positions = new Float32Array(starCount * 3);
        const colors = new Float32Array(starCount * 3);
        
        for (let i = 0; i < starCount * 3; i += 3) {
            // Position
            const radius = 100 + Math.random() * 400;
            const theta = Math.random() * Math.PI * 2;
            const phi = Math.acos((Math.random() * 2) - 1);
            
            positions[i] = radius * Math.sin(phi) * Math.cos(theta);
            positions[i + 1] = radius * Math.sin(phi) * Math.sin(theta);
            positions[i + 2] = radius * Math.cos(phi);
            
            // Color (slight variations)
            const color = new THREE.Color();
            const hue = Math.random() * 0.1 + 0.55; // Blue-white range
            const saturation = Math.random() * 0.2;
            const lightness = 0.6 + Math.random() * 0.4;
            color.setHSL(hue, saturation, lightness);
            
            colors[i] = color.r;
            colors[i + 1] = color.g;
            colors[i + 2] = color.b;
        }
        
        starsGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        starsGeometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        
        const starsMaterial = new THREE.PointsMaterial({
            size: 0.5,
            vertexColors: true,
            transparent: true,
            opacity: 0.8,
            blending: THREE.AdditiveBlending
        });
        
        this.starfield = new THREE.Points(starsGeometry, starsMaterial);
        this.scene.add(this.starfield);
    }
    
    setupPostProcessing() {
        this.composer = new EffectComposer(this.renderer);
        
        // Render pass
        const renderPass = new RenderPass(this.scene, this.camera);
        this.composer.addPass(renderPass);
        
        // Bloom pass for glow effects
        const bloomPass = new UnrealBloomPass(
            new THREE.Vector2(window.innerWidth, window.innerHeight),
            0.5, // strength
            0.4, // radius
            0.85  // threshold
        );
        this.composer.addPass(bloomPass);
    }
    
    async createSaturn() {
        // Create Saturn sphere
        const saturnGeometry = new THREE.SphereGeometry(10, 64, 64);
        
        // Create Saturn texture programmatically
        const canvas = document.createElement('canvas');
        canvas.width = 1024;
        canvas.height = 512;
        const ctx = canvas.getContext('2d');
        
        // Create banded gradient
        const gradient = ctx.createLinearGradient(0, 0, 0, 512);
        gradient.addColorStop(0, '#FAD5A5');
        gradient.addColorStop(0.2, '#E8C39E');
        gradient.addColorStop(0.3, '#D4B896');
        gradient.addColorStop(0.5, '#C9B090');
        gradient.addColorStop(0.7, '#BFA980');
        gradient.addColorStop(0.9, '#B5A070');
        gradient.addColorStop(1, '#AA9560');
        
        ctx.fillStyle = gradient;
        ctx.fillRect(0, 0, 1024, 512);
        
        // Add some storm details
        for (let i = 0; i < 50; i++) {
            const x = Math.random() * 1024;
            const y = Math.random() * 512;
            const size = Math.random() * 20 + 5;
            const opacity = Math.random() * 0.3;
            
            ctx.fillStyle = `rgba(180, 140, 100, ${opacity})`;
            ctx.beginPath();
            ctx.ellipse(x, y, size * 2, size, Math.random() * Math.PI, 0, Math.PI * 2);
            ctx.fill();
        }
        
        const texture = new THREE.CanvasTexture(canvas);
        
        const saturnMaterial = new THREE.MeshPhongMaterial({
            map: texture,
            emissive: 0x332211,
            emissiveIntensity: 0.05,
            shininess: 10
        });
        
        this.saturn = new THREE.Mesh(saturnGeometry, saturnMaterial);
        this.saturn.castShadow = true;
        this.saturn.receiveShadow = true;
        this.saturn.rotation.z = THREE.MathUtils.degToRad(26.73); // Saturn's axial tilt
        this.scene.add(this.saturn);
        
        // Add atmosphere
        const atmosphereGeometry = new THREE.SphereGeometry(10.5, 64, 64);
        const atmosphereMaterial = new THREE.MeshPhongMaterial({
            color: 0xFFE4B5,
            transparent: true,
            opacity: 0.1,
            emissive: 0xFFD700,
            emissiveIntensity: 0.05
        });
        
        this.atmosphere = new THREE.Mesh(atmosphereGeometry, atmosphereMaterial);
        this.scene.add(this.atmosphere);
    }
    
    async createRings() {
        const rings = [
            { inner: 11, outer: 12, color: 0xEEDDCC, opacity: 0.9 },  // D Ring
            { inner: 12.2, outer: 14.5, color: 0xDDCCBB, opacity: 0.8 }, // C Ring
            { inner: 14.7, outer: 18, color: 0xCCBBAA, opacity: 0.7 },   // B Ring
            { inner: 18.5, outer: 20, color: 0xBBAA99, opacity: 0.6 },   // A Ring
            { inner: 20.5, outer: 21, color: 0xAA9988, opacity: 0.4 }    // F Ring
        ];
        
        this.ringSystem = new THREE.Group();
        
        rings.forEach(ring => {
            const geometry = new THREE.RingGeometry(ring.inner, ring.outer, 128, 1);
            const material = new THREE.MeshBasicMaterial({
                color: ring.color,
                side: THREE.DoubleSide,
                transparent: true,
                opacity: ring.opacity
            });
            
            const mesh = new THREE.Mesh(geometry, material);
            mesh.rotation.x = Math.PI / 2;
            mesh.rotation.z = THREE.MathUtils.degToRad(26.73);
            mesh.castShadow = true;
            mesh.receiveShadow = true;
            
            this.ringSystem.add(mesh);
        });
        
        // Add ring particles for detail
        this.createRingParticles();
        
        this.scene.add(this.ringSystem);
    }
    
    createRingParticles() {
        const particleCount = 20000;
        const positions = new Float32Array(particleCount * 3);
        const colors = new Float32Array(particleCount * 3);
        const sizes = new Float32Array(particleCount);
        
        for (let i = 0; i < particleCount; i++) {
            const radius = 11 + Math.random() * 10;
            const angle = Math.random() * Math.PI * 2;
            
            positions[i * 3] = Math.cos(angle) * radius;
            positions[i * 3 + 1] = (Math.random() - 0.5) * 0.1;
            positions[i * 3 + 2] = Math.sin(angle) * radius;
            
            const color = new THREE.Color();
            color.setHSL(0.1, 0.1, 0.7 + Math.random() * 0.3);
            colors[i * 3] = color.r;
            colors[i * 3 + 1] = color.g;
            colors[i * 3 + 2] = color.b;
            
            sizes[i] = Math.random() * 0.5 + 0.1;
        }
        
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        geometry.setAttribute('size', new THREE.BufferAttribute(sizes, 1));
        
        const material = new THREE.PointsMaterial({
            size: 0.1,
            vertexColors: true,
            transparent: true,
            opacity: 0.6,
            blending: THREE.AdditiveBlending
        });
        
        const particles = new THREE.Points(geometry, material);
        particles.rotation.x = Math.PI / 2;
        particles.rotation.z = THREE.MathUtils.degToRad(26.73);
        this.ringSystem.add(particles);
    }
    
    async createMoons() {
        const moonData = [
            { name: 'Mimas', radius: 0.3, distance: 23, period: 0.94, color: 0xE0E0E0 },
            { name: 'Enceladus', radius: 0.4, distance: 26, period: 1.37, color: 0xF8F8FF },
            { name: 'Tethys', radius: 0.8, distance: 29, period: 1.89, color: 0xFFFAFA },
            { name: 'Dione', radius: 0.9, distance: 33, period: 2.74, color: 0xF5F5DC },
            { name: 'Rhea', radius: 1.2, distance: 38, period: 4.52, color: 0xE0E0E0 },
            { name: 'Titan', radius: 2.0, distance: 45, period: 15.95, color: 0xDDB892 },
            { name: 'Iapetus', radius: 1.1, distance: 55, period: 79.32, color: 0x8B7355 }
        ];
        
        this.moonGroup = new THREE.Group();
        
        moonData.forEach(data => {
            const moonGeometry = new THREE.SphereGeometry(data.radius, 32, 32);
            const moonMaterial = new THREE.MeshPhongMaterial({
                color: data.color,
                emissive: data.color,
                emissiveIntensity: 0.02
            });
            
            const moon = new THREE.Mesh(moonGeometry, moonMaterial);
            moon.castShadow = true;
            moon.receiveShadow = true;
            moon.userData = {
                name: data.name,
                distance: data.distance,
                period: data.period,
                angle: Math.random() * Math.PI * 2
            };
            
            // Position moon
            moon.position.x = Math.cos(moon.userData.angle) * data.distance;
            moon.position.z = Math.sin(moon.userData.angle) * data.distance;
            
            this.moons.push(moon);
            this.moonGroup.add(moon);
            
            // Create orbit line
            const orbitPoints = [];
            const segments = 128;
            for (let i = 0; i <= segments; i++) {
                const angle = (i / segments) * Math.PI * 2;
                orbitPoints.push(new THREE.Vector3(
                    Math.cos(angle) * data.distance,
                    0,
                    Math.sin(angle) * data.distance
                ));
            }
            
            const orbitGeometry = new THREE.BufferGeometry().setFromPoints(orbitPoints);
            const orbitMaterial = new THREE.LineBasicMaterial({
                color: 0x444444,
                transparent: true,
                opacity: 0.3
            });
            
            const orbitLine = new THREE.Line(orbitGeometry, orbitMaterial);
            orbitLine.visible = false;
            moon.userData.orbitLine = orbitLine;
            this.moonGroup.add(orbitLine);
        });
        
        this.scene.add(this.moonGroup);
    }
    
    initControls() {
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;
        this.controls.minDistance = 15;
        this.controls.maxDistance = 200;
        this.controls.maxPolarAngle = Math.PI * 0.9;
        this.controls.enablePan = true;
        this.controls.panSpeed = 0.5;
        this.controls.rotateSpeed = 0.5;
    }
    
    initUI() {
        // Play/Pause button
        const playPauseBtn = document.getElementById('play-pause-btn');
        if (playPauseBtn) {
            playPauseBtn.addEventListener('click', () => this.togglePause());
        }
        
        // Reset button
        const resetBtn = document.getElementById('reset-btn');
        if (resetBtn) {
            resetBtn.addEventListener('click', () => this.reset());
        }
        
        // Time speed slider
        const timeSpeedSlider = document.getElementById('time-speed');
        const timeSpeedValue = document.querySelector('.slider-value');
        if (timeSpeedSlider) {
            timeSpeedSlider.addEventListener('input', (e) => {
                this.timeScale = parseFloat(e.target.value);
                if (timeSpeedValue) {
                    timeSpeedValue.textContent = this.timeScale.toFixed(1) + 'x';
                }
            });
        }
        
        // Display toggles
        const toggles = {
            'show-rings': (checked) => {
                this.settings.showRings = checked;
                if (this.ringSystem) this.ringSystem.visible = checked;
            },
            'show-moons': (checked) => {
                this.settings.showMoons = checked;
                if (this.moonGroup) this.moonGroup.visible = checked;
            },
            'show-atmosphere': (checked) => {
                this.settings.showAtmosphere = checked;
                if (this.atmosphere) this.atmosphere.visible = checked;
            },
            'show-orbits': (checked) => {
                this.settings.showOrbits = checked;
                this.moons.forEach(moon => {
                    if (moon.userData.orbitLine) {
                        moon.userData.orbitLine.visible = checked;
                    }
                });
            },
            'show-stats': (checked) => {
                this.settings.showStats = checked;
                if (this.stats) {
                    this.stats.dom.style.display = checked ? 'block' : 'none';
                }
            }
        };
        
        Object.keys(toggles).forEach(id => {
            const checkbox = document.getElementById(id);
            if (checkbox) {
                checkbox.addEventListener('change', (e) => {
                    toggles[id](e.target.checked);
                });
            }
        });
        
        // Camera view buttons
        document.querySelectorAll('.view-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                this.setCameraView(e.target.dataset.view);
            });
        });
        
        // Quality selector
        const qualitySelect = document.getElementById('quality-select');
        if (qualitySelect) {
            qualitySelect.addEventListener('change', (e) => {
                this.setQuality(e.target.value);
            });
        }
        
        // Fullscreen button
        const fullscreenBtn = document.getElementById('fullscreen-btn');
        if (fullscreenBtn) {
            fullscreenBtn.addEventListener('click', () => this.toggleFullscreen());
        }
        
        // Mobile menu toggle
        const mobileMenuToggle = document.getElementById('mobile-menu-toggle');
        if (mobileMenuToggle) {
            mobileMenuToggle.addEventListener('click', () => this.toggleMobileMenu());
        }
    }
    
    setCameraView(view) {
        const views = {
            overview: { position: [30, 15, 30], target: [0, 0, 0] },
            top: { position: [0, 50, 1], target: [0, 0, 0] },
            side: { position: [50, 0, 0], target: [0, 0, 0] },
            rings: { position: [25, 10, 25], target: [0, 0, 0] },
            titan: { position: [45, 5, 45], target: [45, 0, 0] },
            enceladus: { position: [26, 3, 26], target: [26, 0, 0] }
        };
        
        const targetView = views[view];
        if (targetView) {
            // Animate camera to new position
            this.animateCamera(targetView.position, targetView.target);
        }
    }
    
    animateCamera(position, target) {
        // Simple camera animation (could be improved with GSAP or Tween.js)
        const startPos = this.camera.position.clone();
        const endPos = new THREE.Vector3(...position);
        const startTime = Date.now();
        const duration = 1000;
        
        const animate = () => {
            const elapsed = Date.now() - startTime;
            const progress = Math.min(elapsed / duration, 1);
            const eased = this.easeInOutCubic(progress);
            
            this.camera.position.lerpVectors(startPos, endPos, eased);
            this.controls.target.set(...target);
            this.controls.update();
            
            if (progress < 1) {
                requestAnimationFrame(animate);
            }
        };
        
        animate();
    }
    
    easeInOutCubic(t) {
        return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
    }
    
    setQuality(quality) {
        this.settings.quality = quality;
        
        const pixelRatios = {
            low: 1,
            medium: Math.min(window.devicePixelRatio, 1.5),
            high: Math.min(window.devicePixelRatio, 2),
            ultra: window.devicePixelRatio
        };
        
        this.renderer.setPixelRatio(pixelRatios[quality]);
        this.composer.setPixelRatio(pixelRatios[quality]);
        this.onWindowResize();
    }
    
    toggleFullscreen() {
        if (!document.fullscreenElement) {
            document.documentElement.requestFullscreen();
        } else {
            document.exitFullscreen();
        }
    }
    
    toggleMobileMenu() {
        // Implementation for mobile menu
        console.log('Toggle mobile menu');
    }
    
    togglePause() {
        this.isPaused = !this.isPaused;
        const btn = document.getElementById('play-pause-btn');
        if (btn) {
            const icon = btn.querySelector('.btn-icon');
            const text = btn.querySelector('.btn-text');
            if (icon) icon.textContent = this.isPaused ? '▶️' : '⏸️';
            if (text) text.textContent = this.isPaused ? 'Play' : 'Pause';
        }
    }
    
    reset() {
        this.elapsedTime = 0;
        this.timeScale = 1.0;
        
        // Reset moon positions
        this.moons.forEach(moon => {
            moon.userData.angle = Math.random() * Math.PI * 2;
            moon.position.x = Math.cos(moon.userData.angle) * moon.userData.distance;
            moon.position.z = Math.sin(moon.userData.angle) * moon.userData.distance;
        });
        
        // Reset camera
        this.setCameraView('overview');
        
        // Update UI
        const timeSpeedSlider = document.getElementById('time-speed');
        if (timeSpeedSlider) timeSpeedSlider.value = 1;
        const timeSpeedValue = document.querySelector('.slider-value');
        if (timeSpeedValue) timeSpeedValue.textContent = '1.0x';
    }
    
    updateLoadingStatus(status) {
        if (this.loadingStatus) {
            this.loadingStatus.textContent = status;
        }
    }
    
    updateLoadingProgress(percent) {
        if (this.loadingProgress) {
            this.loadingProgress.style.width = percent + '%';
        }
    }
    
    hideLoadingScreen() {
        if (this.loadingScreen) {
            this.loadingScreen.classList.add('fade-out');
            setTimeout(() => {
                this.loadingScreen.style.display = 'none';
            }, 500);
        }
    }
    
    onWindowResize() {
        const width = window.innerWidth;
        const height = window.innerHeight;
        
        this.camera.aspect = width / height;
        this.camera.updateProjectionMatrix();
        
        this.renderer.setSize(width, height);
        this.composer.setSize(width, height);
    }
    
    start() {
        this.isRunning = true;
        this.lastTime = performance.now();
        this.animate();
    }
    
    animate() {
        if (!this.isRunning) return;
        
        requestAnimationFrame(() => this.animate());
        
        const currentTime = performance.now();
        this.deltaTime = (currentTime - this.lastTime) / 1000;
        this.lastTime = currentTime;
        
        if (this.stats) this.stats.begin();
        
        this.update();
        this.render();
        
        if (this.stats) this.stats.end();
    }
    
    update() {
        if (!this.isPaused) {
            const scaledDelta = this.deltaTime * this.timeScale;
            this.elapsedTime += scaledDelta;
            
            // Rotate Saturn
            if (this.saturn) {
                this.saturn.rotation.y += scaledDelta * 0.1;
            }
            
            // Rotate atmosphere slightly slower
            if (this.atmosphere) {
                this.atmosphere.rotation.y += scaledDelta * 0.08;
            }
            
            // Rotate rings
            if (this.ringSystem) {
                this.ringSystem.rotation.z += scaledDelta * 0.05;
            }
            
            // Update moons
            this.moons.forEach(moon => {
                const speed = (1 / moon.userData.period) * scaledDelta;
                moon.userData.angle += speed;
                
                moon.position.x = Math.cos(moon.userData.angle) * moon.userData.distance;
                moon.position.z = Math.sin(moon.userData.angle) * moon.userData.distance;
                
                // Moon rotation
                moon.rotation.y += scaledDelta * 0.5;
            });
            
            // Slowly rotate starfield
            if (this.starfield) {
                this.starfield.rotation.y += scaledDelta * 0.001;
            }
            
            // Update UI stats
            this.updateUIStats();
        }
        
        // Update controls
        if (this.controls) {
            this.controls.update();
        }
    }
    
    updateUIStats() {
        const simTime = document.getElementById('sim-time');
        if (simTime) {
            simTime.textContent = (this.elapsedTime / 86400).toFixed(2) + ' days';
        }
        
        const saturnDay = document.getElementById('saturn-day');
        if (saturnDay) {
            saturnDay.textContent = Math.floor(this.elapsedTime / 38520); // Saturn day in seconds
        }
        
        const activeMoons = document.getElementById('active-moons');
        if (activeMoons) {
            const visibleCount = this.moons.filter(m => m.visible).length;
            activeMoons.textContent = `${visibleCount}/83`;
        }
        
        // Update FPS if performance monitor is visible
        if (this.settings.showStats) {
            const fps = document.getElementById('fps-value');
            if (fps) {
                fps.textContent = Math.round(1 / this.deltaTime);
            }
        }
    }
    
    render() {
        if (this.composer) {
            this.composer.render();
        } else {
            this.renderer.render(this.scene, this.camera);
        }
    }
}

// Initialize simulation when DOM is ready
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', () => {
        new SaturnSimulation();
    });
} else {
    new SaturnSimulation();
}
