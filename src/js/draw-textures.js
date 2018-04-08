import {MeshLambertMaterial} from 'three';


// COLORS
const GREEN = 0x77dd77,

  GREY = 0x888d94,
  CONCRETE = 0xd4d7cd,
  WHITE = 0xfaf8f8,
  FLESH = 0xdfa487,
  BROWN_DARK = 0x472f31,

  BLACK = 0x20261e,
  LIME = 0xabad3d,
  BROWN = 0x41554f,

  YELLOW = 0xf2be54;



export function flatMaterial() {
  return new MeshLambertMaterial( { color: GREEN, flatShading: true} );
}
