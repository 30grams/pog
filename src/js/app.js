import '../css/style.css';
import setup from './game.js';

if (process.env.NODE_ENV !== 'production') {
  console.log('Looks like we are in development mode!');
}

var container = document.createElement('div');
container.id = 'container';

document.body.appendChild(container);

setup();

/* CLEANER LOaDing


async function getComponent() {
  // async import using a Babel pre-processor and the Syntax Dynamic Import Babel Plugin
  var container = document.createElement('div');
  container.id = 'container';

  const THREE = await import(/* webpackChunkName: "threejs" / 'three');
  // dynamic importing to separate a chunk
    
  return container;
}

getComponent().then(component => {
  let element = component(); 
  // Store the element to re-render on hot changes
  document.body.appendChild(element);

  if (module.hot) {
    module.hot.accept('./game.js', function() {
      console.log('Accepting the updated game module!');
      document.body.removeChild(element);
      element = component(); 
      // Re-render the "component" to update the click handler
      document.body.appendChild(element);
    })
  }
})

*/
