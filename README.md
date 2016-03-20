#### pHash
 感知哈希算法的node实现，该算法主要用于提取媒体文件的特征值，可用于识别媒体文件。
但特征值与媒体文件之间并不是唯一对应。一般使用该算法来找出重复的媒体文件。
 
目前这个模块只用于处理图片，使用之前请确保安装了

模块的C++部分使用了[pHash](http://phash.org/)提供的C++库，并在其上做了精简，去除了用不到的部分。

##### Installation
首先请安装[imageMagick](http://www.imagemagick.org/) 或者[GraphicsMagick](http://www.graphicsmagick.org/),
可以使用另一个模块[gm](https://github.com/aheckmann/gm)来验证是否安装正确。

```
  npm install image_phash
```

#### Useage
```javascript

   var image_phash=require("image_phash");
   
   var hash1=image_phash.imageHashSync("./test1.png"),//return the DCT Image Hash
       hash2=image_phash.imageHashSync("./test2.jpg");//maybe you need to storage it for search
	   
   var hamming_distance=image_phash.Hamming_distance(hash1,hash2);//maybe you want use a different threshold to judge the images
   
   var is_similar=image_phash.isSimilar(hash1,hash2);//Threshold set to 26.00.
   
```