// production specificities, merges with common config
const webpack = require('webpack');
const merge = require('webpack-merge');
const common = require('./webpack.common.js');
const UglifyJSPlugin = require('uglifyjs-webpack-plugin');

module.exports = merge(common, {

  mode: 'production',
  devtool: 'source-map',

  plugins: [
    // additional functionalities from webpack
    new UglifyJSPlugin({
      sourceMap: true
    }),
    new webpack.DefinePlugin({
      // system environment variable that Node.js exposes into running scripts
      'process.env.NODE_ENV': JSON.stringify('production')
    })
  ]

});
