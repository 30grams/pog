import * as THREE from 'three';
import Stats from './scripts/stats.min.js';
import populate from './draw-map.js'
import './scripts/TrackballControls.js';

let camera, controls, scene, renderer, stats;
let projector, mouse = {
    x: 0,
    y: 0
  },
  INTERSECTED;
const scale = 1;

export default function setup()
{
  // set up all the 3D objects in the scene   
  init();

  // and let's get cracking!
  animate();
}

// set up all the 3D objects in the scene   
export function init()
{ 
  camera = new THREE.PerspectiveCamera( 30, window.innerWidth / window.innerHeight, 1, 5000 );
  camera.position.set( 0, .3, 7 );

  // WORLD
  scene = new THREE.Scene();
  scene.background = new THREE.Color( 0x000000 );

  // GROUND
  populate(scene);

  let groundGeo = new THREE.PlaneBufferGeometry( 10000, 10000 );
  let groundMat = new THREE.MeshPhongMaterial( { color: 0xffffff, specular: 0x050505 } );
  groundMat.color.setHSL( 0.095, 1, 0.75 );

  let ground = new THREE.Mesh( groundGeo, groundMat );
  ground.rotation.x = -Math.PI/2;
  ground.position.y = -5;
  ground.receiveShadow = true;

  scene.add( ground );

  // LIGHTS
  let hemiLight = new THREE.HemisphereLight( 0xffffff, 0xffffff, 0.6 );
  hemiLight.color.setHSL( 0.6, 1, 0.6 );
  hemiLight.groundColor.setHSL( 0.095, 1, 0.75 );
  hemiLight.position.set( 0, 50, 0 );

  scene.add( hemiLight );

  let hemiLightHelper = new THREE.HemisphereLightHelper( hemiLight, 10 );

  scene.add( hemiLightHelper );


  let dirLight = new THREE.DirectionalLight( 0xffffff, 1 );
  dirLight.color.setHSL( 0.1, 1, 0.95 );
  dirLight.position.set( -1, 1.75, 1 );
  dirLight.position.multiplyScalar( 30 );

  dirLight.castShadow = true;
  dirLight.shadow.mapSize.width = 2048;
  dirLight.shadow.mapSize.height = 2048;

  scene.add( dirLight );

  let d = 10;
  dirLight.shadow.camera.left = -d;
  dirLight.shadow.camera.right = d;
  dirLight.shadow.camera.top = d;
  dirLight.shadow.camera.bottom = -d;

  dirLight.shadow.camera.far = 3500;
  //dirLight.shadow.bias = -0.0001;

  let dirLightHelper = new THREE.DirectionalLightHelper( dirLight, 10 )

  scene.add( dirLightHelper );

  // RENDERER
  renderer = new THREE.WebGLRenderer();
  renderer.setPixelRatio( window.devicePixelRatio );
  renderer.setSize( window.innerWidth, window.innerHeight );

  // SHADOWS
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  // to antialias the shadow

  let container = document.getElementById( 'container' );
  container.appendChild( renderer.domElement );

  // FPS
  stats = new Stats();
  container.appendChild( stats.dom );

  // CONTROLS
  controls = new THREE.TrackballControls( camera );

  controls.rotateSpeed = 1.0;
  controls.zoomSpeed = 1.2;
  controls.panSpeed = 0.8;

  controls.noZoom = false;
  controls.noPan = false;

  controls.staticMoving = true;
  controls.dynamicDampingFactor = 0.3;

  controls.keys = [ 65, 83, 68 ];

  controls.addEventListener( 'change', render );

  // resize
  window.addEventListener( 'resize', onWindowResize, false );

  // when the mouse moves, call the given function
  document.addEventListener('mousemove', onDocumentMouseMove, false);

  // draw
  render();
}

function onWindowResize() {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize( window.innerWidth, window.innerHeight );
  controls.handleResize();
  render();
}

function animate() {
  // main loop
  requestAnimationFrame( animate );
  render();
  update();
}

function onDocumentMouseMove(event) {
  // the following line would stop any other event handler from firing
  // (such as the mouse's TrackballControls)
  // event.preventDefault();

  // update the mouse variable
  mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
  mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
}

function update() {
  // MOUSEOVER

  var vector = new THREE.Vector3(mouse.x, mouse.y, 1);
  vector.unproject(camera);
  var ray = new THREE.Raycaster(camera.position, vector.sub(camera.position).normalize());
  // create a Ray with origin at the mouse position
  //   and direction into the scene (camera direction)

  // create an array containing all objects in the scene with which the ray intersects
  var intersects = ray.intersectObjects(scene.children);

  // INTERSECTED = the object in the scene currently closest to the camera 
  //    and intersected by the Ray projected from the mouse position  

  // if there is one (or more) intersections
  if (intersects.length > 0) {
    // if the closest object intersected is not the currently stored intersection object
    if (intersects[0].object != INTERSECTED) {
      // restore previous intersection object (if it exists) to its original color
      if (INTERSECTED)
        INTERSECTED.material.color.setHex(INTERSECTED.currentHex);
      // store reference to closest object as current intersection object
      INTERSECTED = intersects[0].object;
      // store color of closest object (for later restoration)
      INTERSECTED.currentHex = INTERSECTED.material.color.getHex();
      // set a new color for closest object
      INTERSECTED.material.color.setHex(0xffff00);
    }
  } else // there are no intersections
  {
    // restore previous intersection object (if it exists) to its original color
    if (INTERSECTED)
      INTERSECTED.material.color.setHex(INTERSECTED.currentHex);
    // remove previous intersection object reference
    //     by setting current intersection object to "nothing"
    INTERSECTED = null;
  }

  controls.update();
  stats.update();
}

function render() {
  renderer.render( scene, camera );
}

