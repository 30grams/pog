import {MeshLambertMaterial} from 'three';
import {ShaderMaterial} from 'three';

// COLORS
const GREEN = 0x77dd77, // 119 221 119

  GREY = 0x888d94,
  CONCRETE = 0xd4d7cd,
  WHITE = 0xfaf8f8,
  FLESH = 0xdfa487,
  BROWN_DARK = 0x472f31, // 71 47 49

  BLACK = 0x20261e,
  LIME = 0xabad3d,
  BROWN = 0x41554f,

  YELLOW = 0xf2be54;

// VERTEX SHADER.glsl
const vertexShaderSource = `
varying vec2 vUv;
void main() 
{
  vUv = uv;
  vec4 modelViewPosition = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * modelViewPosition;
}
`;

// FRAGMENT SHADER.glsl
const fragmentShaderSource = `
uniform vec3 topColor;
uniform vec3 baseColor;
varying vec2 vUv;
float frequ = 18.0;
void main() {
  vec3 rgbColor = baseColor;
  float triangleWave = abs(vUv.x*frequ - floor(vUv.x*frequ + 0.5));
  if (vUv.y + 0.1*triangleWave > 0.95) rgbColor = topColor;
  gl_FragColor = vec4(rgbColor, 1.0);
}
`;


// SHADER UNIFORMS
const shaderUniforms = {
  topColor: {value: [119/255, 221/255, 119/255]},
  baseColor: {value: [71/255, 47/255, 49/255]},
};

// MATERIAL SHADER
export function flatMaterialshade() {
  return new ShaderMaterial({
    uniforms: shaderUniforms,
    vertexShader: vertexShaderSource,
    fragmentShader: fragmentShaderSource
  });
}

export function flatMaterial() {
  return [ // An array of materials, which will be used on the subsequent faces
        new ShaderMaterial({
          //lights: true, // crashes...
          uniforms: shaderUniforms,
          vertexShader: vertexShaderSource,
          fragmentShader: fragmentShaderSource,
          flatShading: true
        }),
        new MeshLambertMaterial( { color: GREEN, flatShading: true} )
    ];
  ;
}

export function flatMaterialbak() {
  return new MeshLambertMaterial( { color: GREEN, flatShading: true} );
}
