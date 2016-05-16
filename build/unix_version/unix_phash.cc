#include <node.h>
#include <v8.h>
#include <iostream>
#include "./src/pHash.h"


namespace demo {

using v8::Exception;
using v8::FunctionCallbackInfo;
using v8::Isolate;
using v8::Local;
using v8::Number;
using v8::Object;
using v8::String;
using v8::Value;

std::string ulong64ToString(ulong64 &tar){
    char alpha[10]={'0','1','2','3','4','5','6','7','8','9'};
    int index=0;
    string temp;
    while(tar>10){
        index=tar%10;
        temp=alpha[index]+temp;
        tar=ulong64(tar/10);
    }
     index=tar;
     temp=alpha[index]+temp;
    return temp;
}

ulong64 StringToNumber(std::string tar){
    ulong64 temp=0;
    ulong64 k=1;
    for(int i=tar.length()-1;i>-1;i--){
        temp+=k*(tar[i]-48);
        k*=10;
    }
    return temp;
}

void imageHashSync(const FunctionCallbackInfo<Value>& args){
	Isolate* isolate = args.GetIsolate();
	 if (args.Length()!=1) {
    isolate->ThrowException(Exception::TypeError(
        String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return;
  }
  if(!args[0]->IsString()){
	 isolate->ThrowException(Exception::TypeError(
        String::NewFromUtf8(isolate, "Wrong arguments")));
    return;
  }
  ulong64 tempHash=0;
  v8::String::Utf8Value param1(args[0]->ToString());
  std::string file_path = std::string(*param1);
  int is_error=ph_dct_imagehash(file_path.c_str(),tempHash);
  if(is_error==-1){
	 isolate->ThrowException(Exception::TypeError(
        String::NewFromUtf8(isolate, "Unknow error!")));
    return;  
  }
  std::string ph_result=ulong64ToString(tempHash);
  args.GetReturnValue().Set(String::NewFromUtf8(isolate, ph_result.c_str()));
}
void Hamming_distance(const FunctionCallbackInfo<Value>& args){
	Isolate* isolate = args.GetIsolate();
	 if (args.Length()!=2) {
    isolate->ThrowException(Exception::TypeError(
        String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return;
  }
  if(!args[0]->IsString()||!args[1]->IsString()){
	 isolate->ThrowException(Exception::TypeError(
        String::NewFromUtf8(isolate, "Wrong arguments")));
    return;
  }
   v8::String::Utf8Value param1(args[0]->ToString());
   std::string hash1 = std::string(*param1);
   v8::String::Utf8Value param2(args[1]->ToString());
   std::string hash2 = std::string(*param2);
   int dis_t=ph_hamming_distance(StringToNumber(hash1),StringToNumber(hash2));
   Local<Number> num=Number::New(isolate,dis_t);
   args.GetReturnValue().Set(num);
}

void Init(Local<Object> exports) {
  NODE_SET_METHOD(exports, "imageHashSync", imageHashSync);
  NODE_SET_METHOD(exports, "Hamming_distance", Hamming_distance);
}

NODE_MODULE(unix_phash, Init)

}  