
var the_exact_module=require("./adapt");
var util=require("./util");


exports.imageHashSync=function(image){
    if(util.isString(image)){
        if(util.isExists(image)){
            if(util.isImage(image)){
                try{
                    var hash=the_exact_module.imageHashSync(image);
                    return hash;
                }catch(e){
                    var dl="Failed to get the hash value!Maybe you haven't installed imageMagick or graphicsMagick correctlly!"+
                            "you can use gm(search it on the github) to verify ";
                    throw new Error(dl);
                }
            }else{
                throw new Error("The file '"+image+"' is not a image!");
            }
        }else{
            throw new Error("The file '"+image+"' dosen't exist!");
        }
    }else{
        throw new Error("Wrong argument type!");
    }
}

exports.isSimilar=function(hash1,hash2){
    return the_exact_module.Hamming_distance(hash1,hash2)<26?true:false;
}

exports.Hamming_distance=function(hash1,hash2){
    return the_exact_module.Hamming_distance(hash1,hash2);
}