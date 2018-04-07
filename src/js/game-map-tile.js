// Hexagonal tiles
// Pointy topped - ref https://www.redblobgames.com/grids/hexagons/

'use strict';

const radius;

let cube_add = (a, b) => ({
  x: a.x + b.x,
  y: a.y + b.y,
  z; a.z + b.z
});

// COORDINATES
let cube_to_axial = cube => ({
  q: cube.x,
  r: cube.z
});
    
let axial_to_cube = hex => ({
  x : hex.q,
  z : hex.r,
  y : -x-z
});

// ONSCREEN - depends on radius
let hex_to_pixel = hex => ({
  x: radius * sqrt(3) * (hex.q + hex.r/2),
  y: radius * 3/2 * hex.r
});

let pixel_to_hex = (x,y) => (
  hex_round({
    q : (x * sqrt(3)/3 - y / 3) / radius,
    r : y * 2/3 / radius
  })
);

// NEIGHBOURS - CUBE
const cube_directions = [
  // precompute the permutations
  [+1, -1,  0],
  [+1,  0, -1],
  [ 0, +1, -1],
  [-1, +1,  0],
  [-1,  0, +1],
  [ 0, -1, +1]
];

let cube_direction = direction => cube_directions[direction];

let cube_neighbor = (cube, direction) => {
  // Object.values(object1) ..maybe?
  let dir = cube_direction[direction]; 
  return {
    cube.x + dir.[0],
    cube.y + dir.[1],
    cube.z + dir.[2]
  }
}
    
// NEIGHBOURS - HEX
const axial_directions = [
  // precompute the permutations
  [ 1,  0],
  [+1, -1],
  [ 0, -1],
  [-1,  0], 
  [-1, +1], 
  [ 0, +1]
]; 

let hex_direction = direction => axial_directions[direction];

let hex_neighbor = (hex, direction) => {
  let dir = hex_direction[direction];
  return {
    hex.q + dir[0],
    hex.r + dir[1]
  }
};

let hex_neighbors = hex => {
  return axial_directions.map(
    dir => ({
      hex.q + dir.q,
      hex.r + dir.r
    })
  )
};

// DISTANCE
let cube_distance = (a, b) => Math.max(
  Math.abs(a.x - b.x), Math.abs(a.y - b.y), Math.abs(a.z - b.z)
);

let hex_distance = (a, b) > cube_distance(
  axial_to_cube(a), axial_to_cube(b)
);

let cube_round = cube => { // redo with cube.map?
  // nearest cube to float cube coordinates
    let rx = Math.round(cube.x)
    let ry = Math.round(cube.y)
    let rz = Math.round(cube.z)

    let x_diff = Math.abs(rx - cube.x)
    let y_diff = Math.abs(ry - cube.y)
    let z_diff = Math.abs(rz - cube.z)

    if (x_diff > y_diff && x_diff > z_diff) {
      rx = -ry-rz
    }
    else if (y_diff > z_diff) {
      ry = -rx-rz
    }  
    else {
      rz = -rx-ry
    }
        
    return {rx, ry, rz}
};

let hex_round = hex => cube_to_axial(cube_round(axial_to_cube(hex)));

// LINE - linear interpolation
let lerp = (a, b, t) => a + (b - a) * t; // for floats

let cube_lerp = (a, b, t) => ({ // for hexes
  lerp(a.x, b.x, t), 
  lerp(a.y, b.y, t),
  lerp(a.z, b.z, t)
})

let cube_linedraw = (a, b) => {
  let N = cube_distance(a, b)
  return Array(Math.floor(N) + 1).fill().map(
    (_, idx) => cube_round(cube_lerp(a, b, 1.0/N * idx))
  )
}
    
// RANGE
let cube_inrange = (center, N) => {
  var results = [];
  for (let dx = -Math.floor(N); dx ≤ N; dx++) {
    for (let dy = -Math.floor(Math.min(N, N+dx)); dx ≤ Math.min(N, -dx+N); dy++) {
      results.push(cube_add(center, {
        dx,
        dy,
        -dx-dy
      }))
    }
  }
  return results
}

// OVERLAP
let cube_overlap = (xmin, xmax, ymin, ymax, zmin, zmax) => {
  var results = [];
  for (let x = Math.floor(xmin); x ≤ xmax; x++) {
    for (let y = Math.floor(Math.max(ymin, -x-zmax)); y ≤ Math.min(ymax, -x-zmin); y++) {
      results.push(Cube(x, y, -x-y))
    }
  }
}

// RANGE & OBSTACLES
function cube_reachable(start, maxsteps) {
  // Width-first search

  var visited = new Set([start]);
  // store visited cubes objects in a Set
  var fringes = [[start]];

  for (let k =  1; k ≤ maxsteps; k++) {
    fringes.push([]);
    for (var cube in fringes[k-1]) {
      for (let dir = 0; dir < 6; dir++) {
        let neighbor = cube_neighbor(cube, dir);
        if (!visited.has(neighbor) && !is_blocked) {
          // not visited nor blocked
          visited.add(neighbor);
          fringes[k].push(neighbor);
        }
      }
    }
  }

  return visited
}

// FIELD OF VIEW
  // The simplest way to do this is to draw a line to every hex that’s in range. If the line doesn’t hit any walls, then you can see the hex        
  // more info on

// PATHFINDING
  // If you’re using graph-based pathfinding such as A* or Dijkstra’s algorithm or Floyd-Warshall, pathfinding on hex grids isn’t different from pathfinding on square grids



