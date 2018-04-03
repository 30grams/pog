// webpack config file
const path = require('path');
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const CleanWebpackPlugin = require('clean-webpack-plugin');
const ExtractTextPlugin = require("extract-text-webpack-plugin");


module.exports = {

  entry: './src/js/app.js',
  // Here the application starts executing
  // and webpack starts bundling
  output: {
     // options related to how webpack emits results
    path: path.resolve(__dirname, "dist"),
    filename: '[name].bundle.js',
    chunkFilename: '[name].bundle.js'
  },

  plugins: [
    // additional functionalities from webpack
    new ExtractTextPlugin({
      disable: process.env.NODE_ENV !== 'production',
      filename: "styles.css"
    }),
    new CleanWebpackPlugin(['dist']),
    new HtmlWebpackPlugin({
      title: 'PoG'
    })
  ],

  optimization: {
    namedModules: true,
    runtimeChunk: true, 
    // to avoid all hashes of generated file change every time a piece of code change in 1 file ... (and spare 4kb)
    splitChunks: {
      // split the vendor and application code
      chunks: 'all',
      name : false,
      cacheGroups: {
        commons: { 
          test: /[\\/]node_modules[\\/]/,
          name: "vendors", 
          chunks: "all" 
        }
      }
    }
  },

  module:{
    // configuration regarding modules
    rules: [
    // rules for modules (configure loaders, parser options, etc.)
      {
        test: /\.js$/,
        include: [
          path.resolve(__dirname, 'src/js')
        ],
        exclude: [
          path.resolve(__dirname, 'src/js/scripts')
        ],
        // these are matching conditions, each accepting a regular expression or string
        // test and include have the same behavior, both must be matched
        // exclude must not be matched (takes preferrence over test and include)
        loader: 'babel-loader',
        // the loader which should be applied, it'll be resolved relative to the context
        // -loader suffix is no longer optional in webpack2 for clarity reasons
        // see webpack 1 upgrade guide
        options: {
          presets: ['es2015']
        }
        // options for the loader
      },

      {
        test: /\.css$/,
        // loader for css files
        include: [
          path.resolve(__dirname, 'src/css')
        ],
        use: ExtractTextPlugin.extract({
          // separates css in a file rather than inline
          fallback: 'style-loader',
          use: 'css-loader'
        })
      },

      {
        test: /\.(woff|woff2|eot|ttf|otf)$/,
        // loader for the font files
        loader: 'file-loader'
      }

    ]
  }
};
