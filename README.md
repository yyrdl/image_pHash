#### pHash
 感知哈希算法的node实现，该算法主要用于提取媒体文件的特征值，可用于识别媒体文件。
但特征值与媒体文件之间并不是唯一对应。一般使用该算法来找出重复的媒体文件。

模块的C++部分使用了[pHash](http://phash.org/)提供的C++库，并在其上做了精简，去除了用不到的部分。

目前这个模块只用于处理图片

  A perceptual hash is a fingerprint of a multimedia file derived from various features from its content.
You can use it to find simillar image.
  The main C++ code is from [pHash](http://phash.org/),and I deleted some code that will not be used.Such as 
the code of trying to load image by built-in jpeglib ,and the code that deal with text and video.etc. This module
mainly deals with images .
  And there is still a problem,that is the 'CImg.h' try to throw an error when it failed to load the image by external
service ,which is surpported by [imageMagick](http://www.imagemagick.org/) or [GraphicsMagick](http://www.graphicsmagick.org/).
It will kill the process.I have tried to catch the error using c++ code or javascript code, but failed!
The only way to solve this is to modify the "CImg.h" file,but there is more than 50 thousand lines .......

  Anyway, this module works well  if you have installed [imageMagick](http://www.imagemagick.org/) or [GraphicsMagick](http://www.graphicsmagick.org/) correctly.
And I have tested this module in windows(win10)  and linux(ubutun).
  


#### Installation
首先请安装[imageMagick](http://www.imagemagick.org/) 或者[GraphicsMagick](http://www.graphicsmagick.org/),
可以使用另一个模块[gm](https://github.com/aheckmann/gm)来验证是否安装正确。

 At first ,please install [imageMagick](http://www.imagemagick.org/) or [GraphicsMagick](http://www.graphicsmagick.org),you can use [gm](https://github.com/aheckmann/gm) to verify the 
installation.

```
  npm install image_phash
```

#### Useage
```javascript

   var image_phash=require("image_phash");
   
   var hash1=image_phash.imageHashSync("./test1.png"),//return the DCT Image Hash
       hash2=image_phash.imageHashSync("./test2.jpg");//maybe you need to storage 
	   //it for search
	   
   var hamming_distance=image_phash.Hamming_distance(hash1,hash2);//maybe you want 
   //use a different threshold to judge the images
   
   var is_similar=image_phash.isSimilar(hash1,hash2);//Threshold set to 26.00. this 
   //function return a boolean value,while 'true' means these two image is similar! 
```

  
