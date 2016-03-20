#### pHash
 ��֪��ϣ�㷨��nodeʵ�֣����㷨��Ҫ������ȡý���ļ�������ֵ��������ʶ��ý���ļ���
������ֵ��ý���ļ�֮�䲢����Ψһ��Ӧ��һ��ʹ�ø��㷨���ҳ��ظ���ý���ļ���
 
Ŀǰ���ģ��ֻ���ڴ���ͼƬ��ʹ��֮ǰ��ȷ����װ��

ģ���C++����ʹ����[pHash](http://phash.org/)�ṩ��C++�⣬�����������˾���ȥ�����ò����Ĳ��֡�

##### Installation
�����밲װ[imageMagick](http://www.imagemagick.org/) ����[GraphicsMagick](http://www.graphicsmagick.org/),
����ʹ����һ��ģ��[gm](https://github.com/aheckmann/gm)����֤�Ƿ�װ��ȷ��

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