// Hexagonal tiles
// Pointy topped - ref https://www.redblobgames.com/grids/hexagons/

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

// NEIGHBOURS
const axial_directions = [
   {1,  0}, {+1, -1}, { 0, -1},
   {-1, 0}, {-1, +1}, { 0, +1}
];   // precompute the permutations

let hex_direction = direction => axial_directions[direction];

let hex_neighbor = (hex, direction) => {
  let dir = hex_direction[direction];
  return {hex.q + dir.q, hex.r + dir.r}
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

let hex_distance = (a,b) > cube_distance(
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
  var results = [] // range
  for (let dx = -Math.floor(N), dx ≤ N, dx++) {
    for (let dy = -Math.floor(Math.min(N, N+dx)), dx ≤ Math.min(N, -dx+N), dy++) {
      results.append(cube_add(center, {
        dx,
        dy,
        -dx-dy
      }))
    }
  }
  return results
}

var results = [] // overlap
for each xmin ≤ x ≤ xmax:
    for each max(ymin, -x-zmax) ≤ y ≤ min(ymax, -x-zmin):
        var z = -x-y
        results.append(Cube(x, y, z))

function cube_reachable(start, movement): // range with obstacles
    var visited = set()
    add start to visited
    var fringes = []
    fringes.append([start])

    for each 1 < k ≤ movement:
        fringes.append([])
        for each cube in fringes[k-1]:
            for each 0 ≤ dir < 6:
                var neighbor = cube_neighbor(cube, dir)
                if neighbor not in visited, not blocked:
                    add neighbor to visited
                    fringes[k].append(neighbor)

    return visited

// FIELD OF VIEW
  // The simplest way to do this is to draw a line to every hex that’s in range. If the line doesn’t hit any walls, then you can see the hex        
  // more info on

// PATHFINDING
  // If you’re using graph-based pathfinding such as A* or Dijkstra’s algorithm or Floyd-Warshall, pathfinding on hex grids isn’t different from pathfinding on square grids

// Map storage


