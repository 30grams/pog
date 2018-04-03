// developpment specificities, merges with common config
const merge = require('webpack-merge');
const common = require('./webpack.common.js');
const webpack = require('webpack');

module.exports = merge(common, {

  devtool: 'inline-source-map',
  devServer: {
    contentBase: './dist'
    },
  mode: 'development',
  
  plugins: [
    // additional functionalities from webpack
    new webpack.NamedModulesPlugin(),
    new webpack.HotModuleReplacementPlugin(),
  ]

});
