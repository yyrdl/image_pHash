/**
 * Created by yyrdl
 */

var platForm=require("os").platform();
var phash;
if(platForm=="win32"){
  phash=require("../Release/win32/win32_phash.node");
}else{
  phash=require("../Release/unix/unix_phash.node");
}

module.exports=phash;


