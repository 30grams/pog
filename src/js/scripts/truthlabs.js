 var modelList = ['cockpit1.js', 'cockpit2.js', 'fin1.js', 'fin2.js', 'wing1.js', 'wing2.js', 'rocket1.js', 'rocket2.js'];
 var loadCount = 0;
 var loaded = false;
 var path = 'https://s3-us-west-2.amazonaws.com/experiments.truthlabs.com/LowPoly/models/';
 //var path = '/models/';

 function initObjects() {

     group = new THREE.Object3D();
    // group.position.set(0, 0, 150);

     loadShip();
     loadGround();
     loadSky();
     loadPlanet();
     loadStars();
     loadLights();
 }

 function loadShip() {

     var loader = new THREE.JSONLoader();
     loader.load('https://s3-us-west-2.amazonaws.com/experiments.truthlabs.com/LowPoly/models/ship.js', function (geometry, materials) {

         var material = new THREE.MeshLambertMaterial({
             ambient: 0x999999,
             color: 0x777777,
             emissive: 0x685b6e

         });

         ship = new THREE.Mesh(geometry, material);
         //ship.position.set, -5, 0);
         ship.scale.multiplyScalar(10);

         ship.castShadow = true;
         ship.receiveShadow = false;
         group.add(ship);


     });
 }

 function loadStars() {


     var stars = new THREE.CylinderGeometry(1000, 1000, 500, 50, 50, true);
     var star_stuff = new THREE.ParticleBasicMaterial({
         color: 0xffffff,
         size: 6,
         opacity: .66,
         blending: THREE.AdditiveBlending,
         transparent: true
     });

     star_system = new THREE.ParticleSystem(stars, star_stuff);
     var vertices = star_system.geometry.vertices;

     for (i = 0; i < vertices.length; i++) {
         vertices[i].x += Math.floor(Math.random() * 100);
         vertices[i].y += Math.floor(Math.random() * 100);
         vertices[i].z += Math.floor(Math.random() * 500);
     }

     scene.add(star_system);

 }

 function loadGround() {

     var loader = new THREE.JSONLoader();
     loader.load('https://s3-us-west-2.amazonaws.com/experiments.truthlabs.com/LowPoly/models/scene.js', function (geometry, materials) {

         var material = new THREE.MeshLambertMaterial({
             ambient: 0x999999,
             color: 0xfd7a7b,
             emissive: 0x6a2323
         });

         ground = new THREE.Mesh(geometry, material);
         ground.rotation.set(-55, -Math.PI / 1, 4.9);
         ground.position.set(20, -30, -150);
         ground.scale.multiplyScalar(.4);

         ground.castShadow = false;
         ground.receiveShadow = true;
         scene.add(ground);

         loaded = true;
         var loadingDom = document.getElementById('loading');
         loading.style.display = 'none';
         group.scale.multiplyScalar(1.4);
         scene.add(group);

         animate();

     });
 }


 function loadPlanet() {

     var loader = new THREE.JSONLoader();
     loader.load('https://s3-us-west-2.amazonaws.com/experiments.truthlabs.com/LowPoly/models/planet.js', function (geometry, materials) {

         var material = new THREE.MeshLambertMaterial({
             ambient: 0x999999,
             color: 0xfd7a7b,
             emissive: 0x6a2323            
         });

         planet = new THREE.Mesh(geometry, material);
         planet.rotation.set(-55.1, -Math.PI / 1, 9.3);
         planet.position.set(-50, 30, -450);
         planet.scale.multiplyScalar(.3);

         planet.castShadow = false;
         planet.receiveShadow = true;
         scene.add(planet);

     });
 }

 function loadSky() {

     // SKYDOME

     var vertexShader = document.getElementById('vertexShader').textContent;
     var fragmentShader = document.getElementById('fragmentShader').textContent;
     var uniforms = {
         topColor: {
             type: "c",
             value: new THREE.Color(0x562e35)
         },
         bottomColor: {
             type: "c",
             value: new THREE.Color(0x382630)
         },
         offset: {
             type: "f",
             value: 30
         },
         exponent: {
             type: "f",
             value: .5
         }
     }

     var skyGeo = new THREE.SphereGeometry(1000, 32, 15);
     var skyMat = new THREE.ShaderMaterial({
         vertexShader: vertexShader,
         fragmentShader: fragmentShader,
         uniforms: uniforms,
         side: THREE.BackSide
     });
     var material = new THREE.MeshBasicMaterial({
         color: 0xff0000
     });
     var sky = new THREE.Mesh(skyGeo, skyMat);
     sky.position.set(0, 0, 0);
     scene.add(sky);

 }

 function loadLights() {

     // LIGHTS

     var light = new THREE.SpotLight(0x999999, 1);
     light.castShadow = true;
     light.shadowDarkness = .3;
     light.position.set(0, 330, 0);
     scene.add(light);

     var light = new THREE.PointLight(0x999999, .6);
     light.position.set(500, 200, 0);
     scene.add(light);

     hemiLight = new THREE.HemisphereLight(0xff7a81, 0x000000, .5);
     hemiLight.position.y = 500;
     scene.add(hemiLight);


 }