#### pHash ![build status](https://travis-ci.org/yyrdl/image_pHash.svg?branch=master)
 感知哈希算法的node实现，该算法主要用于提取媒体文件的特征值，可用于识别媒体文件。
但特征值与媒体文件之间并不是唯一对应。一般使用该算法来找出重复的媒体文件。

模块的C++部分使用了[pHash](http://phash.org/)提供的C++库，并在其上做了精简，去除了用不到的部分。

目前这个模块只用于处理图片

  A perceptual hash is a fingerprint of a multimedia file derived from various features from its content.
You can use it to find simillar image.
  The main C++ code is from [pHash](http://phash.org/),and I deleted some code that will not be used.Such as the code of trying to load image by built-in jpeglib ,and the code that deal with text and video.etc. .And I removed the error that throw out by 'CImg.h',This change avoids the crash of your process when the CImg Object try to throw an error!
   
  This module mainly deals with images and  works well  if you have installed [imageMagick](http://www.imagemagick.org/) or [GraphicsMagick](http://www.graphicsmagick.org/) correctly.
And I have tested this module in windows(win10)  and linux(ubutun).
  


#### Installation
首先请安装[imageMagick](http://www.imagemagick.org/) 或者[GraphicsMagick](http://www.graphicsmagick.org/),
可以使用另一个模块[gm](https://github.com/aheckmann/gm)来验证是否安装正确。

 At first ,please install [imageMagick](http://www.imagemagick.org/) or [GraphicsMagick](http://www.graphicsmagick.org),you can use [gm](https://github.com/aheckmann/gm) to verify the 
installation.

```
  npm install image_phash
```
the version in the npm warehouse only support node 4.*.For other node version ,please rebuild.
#### Build
如果相应的编译环境设置正确，运行下面的命令即可自动编译C++扩展

if you have set the compiling environment of special os correctly ,just run the command line below ,it will build the C++ addon  automatically.
##### command line
```
   node build.js
```
##### compiling environment
编译环境
###### unix
  * node-gyp 
  * python2.7 
  * g++ 
  
###### windows
  * node-gyp 
  * python2.7
  * windows的C++编译环境，安装visual studio会自动安装相应的编译环境
  
#### Useage
```javascript

   var image_phash=require("image_phash"),
       path=require("path");
   
   var hash1=image_phash.imageHashSync(path.join(__dirname,"./test1.png")),//return the DCT Image Hash
       hash2=image_phash.imageHashSync(path.join(__dirname,"./test2.jpg"));//maybe you need to storage
	   //it for search
	   
   var hamming_distance=image_phash.Hamming_distance(hash1,hash2);//maybe you want 
   //use a different threshold to judge the images
   
   var is_similar=image_phash.isSimilar(hash1,hash2);//Threshold set to 26.00. this 
   //function return a boolean value,while 'true' means these two image is similar! 
```

  
