// by yyrdl 
var fs=require("fs");
var img_extname_list = {
	"bmp" : true,
	"jpg" : true,
	"jpeg" : true,
	"jpe" : true,
	"jfif" : true,
	"jif" : true,
	"png" : true,
	"ppm" : true,
	"pgm" : true,
	"pnm" : true,
	"pbm" : true,
	"pnk" : true,
	"pfm" : true,
	"tif" : true,
	"tiff" : true,
	"exr" : true,
	"cr2" : true,
	"crw" : true,
	"dcr" : true,
	"mrw" : true,
	"nef" : true,
	"orf" : true,
	"pix" : true,
	"ptx" : true,
	"raf" : true,
	"srf" : true,
	"gif" : true
};

function ext_name(file_path){
	file_path=file_path.split(" ").join();
	var temp=file_path.split(".");
	return temp[temp.length-1];
}

exports.isImage=function(file_path){
	return img_extname_list[ext_name(file_path)]?true:false;
}
exports.isExists=function(file_path){
	return fs.existsSync(file_path)
}

exports.isString=function(str){
	return Object.prototype.toString.call(str)==="[object String]";
}

exports.isNumberCharacterSequence=function(tar){
	if(!exports.isString(tar)){
		return false;
	}
	return isNaN(parseInt(tar))?false:true;
}
