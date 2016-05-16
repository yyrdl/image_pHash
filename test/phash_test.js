/**
 * Created by jason on 2016/5/16.
 */
var expect=require("chai").expect;
var phash=require("../lib");

describe("Should be simillar with itself",function(){
    var hash1,hash2,hash3,hash4,hash5;
    before(function(){
            hash1=phash.imageHashSync("./test_images/1.jpg"),
            hash2=phash.imageHashSync("./test_images/2.jpg"),
            hash3=phash.imageHashSync("./test_images/3.jpg"),
            hash4=phash.imageHashSync("./test_images/4.jpg"),
            hash5=phash.imageHashSync("./test_images/5.jpg");
    });
    it("1.jpg should be simillar with 1.jpg",function(){

        var result=phash.isSimilar(hash1,hash1);
        expect(result).to.equal(true);
    });

    it("2.jpg should be simillar with 2.jpg",function(){
        var result=phash.isSimilar(hash2,hash2);
        expect(result).to.equal(true);
    });

    it("3.jpg should be simillar with 3.jpg",function(){
        var result=phash.isSimilar(hash3,hash3);
        expect(result).to.equal(true);
    });

    it("4.jpg should be simillar with 4.jpg",function(){
        var result=phash.isSimilar(hash4,hash4);
        expect(result).to.equal(true);
    });


    it("5.jpg should be simillar with 5.jpg",function(){
        var result=phash.isSimilar(hash5,hash5);
        expect(result).to.equal(true);
    })
});

describe("Should be simillar",function(){

    it("1.jpg should be simillar with 2.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/1.jpg"),
            hash2=phash.imageHashSync("./test_images/2.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    })

    it("1.jpg should be simillar with 3.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/1.jpg"),
            hash2=phash.imageHashSync("./test_images/3.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    })

    it("1.jpg should be simillar with 4.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/1.jpg"),
            hash2=phash.imageHashSync("./test_images/4.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    })

    it("2.jpg should be simillar with 3.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/2.jpg"),
            hash2=phash.imageHashSync("./test_images/3.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    })

    it("2.jpg should be simillar with 4.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/2.jpg"),
            hash2=phash.imageHashSync("./test_images/4.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    })

    it("3.jpg should be simillar with 4.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/3.jpg"),
            hash2=phash.imageHashSync("./test_images/4.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    })
});

describe("Should not be simillar!",function(){

    it("1.jpg should not be simillar with 5.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/1.jpg"),
            hash2=phash.imageHashSync("./test_images/5.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    });

    it("2.jpg should be simillar with 5.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/2.jpg"),
            hash2=phash.imageHashSync("./test_images/5.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    });

    it("3.jpg should be simillar with 5.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/3.jpg"),
            hash2=phash.imageHashSync("./test_images/5.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    });

    it("4.jpg should be simillar with 5.jpg",function(){
        var hash1=phash.imageHashSync("./test_images/4.jpg"),
            hash2=phash.imageHashSync("./test_images/5.jpg");
        var result=phash.isSimilar(hash1,hash2);
        expect(result).to.equal(true);
    });
})