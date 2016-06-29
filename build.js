/**
 * Created by yyrdl on 2016/6/29.
 */
"use strict"
const exec=require("child_process").exec;
const path=require("path");
const fs=require("fs");
const os=require("os");

const buildWinVersion=function(){
     let  workdirectory=path.join(__dirname,"./build/win32");

     let config_process=exec("node-gyp configure",{
         cwd:workdirectory
     },function(err,stdout,stderr){
         if(err){
             console.error(err.stack);
             console.log(stderr);
         }
     });
    let transform_node_addon=function(){
       let addon_path=path.join(__dirname,"./lib/Release/win32");
        if(fs.existsSync(addon_path+"/win32_phash.node")){
            fs.unlinkSync(addon_path+"/win32_phash.node");
        }
       if(fs.existsSync(path.join(workdirectory,"./build/Release/win32_phash.node"))){
           fs.linkSync(path.join(workdirectory,"./build/Release/win32_phash.node"),path.join(addon_path,"./win32_phash.node"));
           console.log("Build Done!");

       }else{
           throw new Error("Build Failed!");
       }
    };

    let nextStep=function(){
       let build_process=exec("node-gyp build",{
           cwd:workdirectory
       },function(err,stdout,stderr){
           if(err){
               console.error(err.stack);
               console.log(stderr);
           }
       });
        build_process.on("error",function(err){
            throw err;
        });

        build_process.on("exit",function(){

               transform_node_addon();
        });
    };
    config_process.on("error",function(err){
        throw err;
    })
    config_process.on("exit",function(){
            nextStep();
    });
}

const buildUnixVersion=function(){
    let  workdirectory=path.join(__dirname,"./build/unix_version");

    let config_process=exec("node-gyp configure",{
        cwd:workdirectory
    },function(err,stdout,stderr){
        if(err){
            console.error(err.stack);
            console.log(stderr);
        }
    });
    let transform_node_addon=function(){
        let addon_path=path.join(__dirname,"./lib/Release/unix");
        if(fs.existsSync(addon_path+"/unix_phash.node")){
            fs.unlinkSync(addon_path+"/unix_phash.node");
        }
        if(fs.existsSync(path.join(workdirectory,"./build/Release/unix_phash.node"))){
            fs.linkSync(path.join(workdirectory,"./build/Release/unix_phash.node"),path.join(addon_path,"./unix_phash.node"));
            console.log("Build Done!");
        }else{
            throw new Error("Build Failed!");
        }
    };
    let fexceptions=function(){
        var file=path.join(workdirectory,"./build/unix_phash.target.mk");
        var unix_phash_target_mk=fs.readFileSync(file);
        unix_phash_target_mk=unix_phash_target_mk.toString().replace(/fno-exceptions/g,"fexceptions");
        fs.unlinkSync(file);
        fs.appendFile(file,unix_phash_target_mk);
    }
    let nextStep=function(){
        let build_process=exec("node-gyp build",{
            cwd:workdirectory
        },function(err,stdout,stderr){
            if(err){
                console.error(err.stack);
                console.log(stderr);
            }
        });
        build_process.on("error",function(err){
            throw err;
        });

        build_process.on("exit",function(){

            transform_node_addon();
        });
    };
    config_process.on("error",function(err){
        throw err;
    })
    config_process.on("exit",function(){
        fexceptions();
        nextStep();
    });
}

const build=function(){
    var os_type=os.platform();
    if(os_type=="win32"){
        buildWinVersion();
    }else{
        buildUnixVersion();
    }
}

build();

