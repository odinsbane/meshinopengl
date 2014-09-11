#include "TiffWriter.h"

TiffWriter::TiffWriter(const char* filename, int height, int width){
    
    image = TIFFOpen(filename, "w");
    this->height = height;
    this->width = width;
    count=0;
    bpp=3;
    
}

void TiffWriter::writeFrame(char* buffer){
    // We need to set some values for basic tags before we can add any data
    TIFFSetField(image, TIFFTAG_IMAGEWIDTH, width );
    TIFFSetField(image, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, bpp);

    TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    //TIFFSetField(image, TIFFTAG_ORIENTATION, ORIENTATION_LEFTBOT);
    // Write the information to the file
    
    for( int j = 0; j< height; j++){
        TIFFWriteScanline(image, &buffer[j*width*bpp], j, 0);
    }
    
    TIFFFlush(image);

    count++;
}

void TiffWriter::close(){
    TIFFClose(image);
}

int TiffWriter::getCount(){return count;};