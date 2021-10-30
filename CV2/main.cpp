#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
using namespace std;

#define UnknownError 0
#define ImageNotExists 1
#define BinaryFileError 2
#define UnrecognizedCompressionType 3

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;

inline void errorMessage(uint8 errorCode)
{
    switch (errorCode) {
        case UnknownError:
            puts("Unknown Error!");
            break;
        case ImageNotExists:
            puts("Image not exists!");
            break;
        case BinaryFileError:
            puts("The binary file error!");
            break;
        case UnrecognizedCompressionType:
            puts("The compression type is unrecognized!");
            break;
        default:
            puts("Error code is miss.");
            break;
    }
    exit(0);
}

struct RGB {
    double r, g, b;
};

struct BitMapFileHeader {
    uint16 bfType;           // File Type ("BM" required)
    uint32 bfSize;           // Size of File
    uint16 bfReserved1;      // Reserved ("00 00" required)
    uint16 bfReserved2;      // Reserved ("00 00" required)
    uint32 bfOffbits;        // Offset Bits from File Header to Image Data
    inline void read(FILE *imgin) {
        fread(&bfType, 2, 1, imgin);
        fread(&bfSize, 4, 1, imgin);
        fread(&bfReserved1, 2, 1, imgin);
        fread(&bfReserved2, 2, 1, imgin);
        fread(&bfOffbits, 4, 1, imgin);
    }
    inline void write(FILE *imgout) {
        fwrite(&bfType, 2, 1, imgout);
        fwrite(&bfSize, 4, 1, imgout);
        fwrite(&bfReserved1, 2, 1, imgout);
        fwrite(&bfReserved2, 2, 1, imgout);
        fwrite(&bfOffbits, 4, 1, imgout);
    }
};

struct BitMapInfoHeader {
    uint32 biSize;           // Size of Info Header
    uint32 biWidth;          // Width of Image (px)
    uint32 biHeight;         // Height of Image (px)
    uint16 biPlanes;         // Number of Planes with Color ("01 00" required)
    uint16 biBitCount;       // Bits per Pixel
    uint32 biCompression;    // Compression Type ("00 00 00 00": None, "01 00 00 00": RLE8)
    uint32 biSizeImage;      // Size of Image (can be "00 00 00 00" if biCompression is "00 00 00 00")
    uint32 biXPelsPerMeter;  // Horizontal Resolution (pixel per meter, same as above)
    uint32 biYPelsPerMeter;  // Vertical Resolution (pixel per meter, same as above)
    uint32 biClrUsed;        // Number of Used Indexes
    uint32 biClrImportant;   // Number of Important Indexes
    inline void read(FILE *imgin) {
        fread(&biSize, 4, 1, imgin);
        fread(&biWidth, 4, 1, imgin);
        fread(&biHeight, 4, 1, imgin);
        fread(&biPlanes, 2, 1, imgin);
        fread(&biBitCount, 2, 1, imgin);
        fread(&biCompression, 4, 1, imgin);
        fread(&biSizeImage, 4, 1, imgin);
        fread(&biXPelsPerMeter, 4, 1, imgin);
        fread(&biYPelsPerMeter, 4, 1, imgin);
        fread(&biClrUsed, 4, 1, imgin);
        fread(&biClrImportant, 4, 1, imgin);
    }
    inline void write(FILE *imgout) {
        fwrite(&biSize, 4, 1, imgout);
        fwrite(&biWidth, 4, 1, imgout);
        fwrite(&biHeight, 4, 1, imgout);
        fwrite(&biPlanes, 2, 1, imgout);
        fwrite(&biBitCount, 2, 1, imgout);
        fwrite(&biCompression, 4, 1, imgout);
        fwrite(&biSizeImage, 4, 1, imgout);
        fwrite(&biXPelsPerMeter, 4, 1, imgout);
        fwrite(&biYPelsPerMeter, 4, 1, imgout);
        fwrite(&biClrUsed, 4, 1, imgout);
        fwrite(&biClrImportant, 4, 1, imgout);
    }
};

struct BitMapColorIndex {
    uint8 r, g, b, alpha;
    inline void read(FILE *imgin) {
        fread(&b, 1, 1, imgin);
        fread(&g, 1, 1, imgin);
        fread(&r, 1, 1, imgin);
    }
    inline void write(FILE *imgout) {
        fwrite(&b, 1, 1, imgout);
        fwrite(&g, 1, 1, imgout);
        fwrite(&r, 1, 1, imgout);
    }
};

class Image {
private:
    BitMapFileHeader fileHeader;
    BitMapInfoHeader infoHeader;
    BitMapColorIndex *index;
    int **pixelIndex;
public:
    void open(char *s) {
        FILE *imgin = fopen(s, "rb");
        if (imgin == 0) errorMessage(ImageNotExists);
        fileHeader.read(imgin);
        infoHeader.read(imgin);
        index = (BitMapColorIndex*)malloc(sizeof(BitMapColorIndex) * infoHeader.biHeight * infoHeader.biWidth);
        pixelIndex = (int**)malloc(sizeof(int*) * infoHeader.biHeight);
        int tot = 0;
        for (int i = infoHeader.biHeight - 1; i >= 0; --i) {
            pixelIndex[i] = (int*)malloc(sizeof(int) * infoHeader.biWidth);
            for (int j = 0; j < infoHeader.biWidth; ++j) {
                index[tot].read(imgin);
                pixelIndex[i][j] = tot++;
            }
            int lineByte = (infoHeader.biWidth * infoHeader.biBitCount / 8 + 3) / 4 * 4;
            for (int j = infoHeader.biWidth * 3; j < lineByte; ++j) {
                uint8 tmp;
                fread(&tmp, 1, 1, imgin);
            }
        }
        fclose(imgin);
        return;
    }
    void save(char *s) {
        FILE *imgout = fopen(s, "wb");
        fileHeader.write(imgout);
        infoHeader.write(imgout);
        uint8 empty = 0;
        for (int i = infoHeader.biHeight - 1; i >= 0; --i) {
            for (int j = 0; j < infoHeader.biWidth; ++j) {
                index[pixelIndex[i][j]].write(imgout);
            }
            int lineByte = (infoHeader.biWidth * infoHeader.biBitCount / 8 + 3) / 4 * 4;
            for (int j = infoHeader.biWidth * 3; j < lineByte; ++j) {
                fwrite(&empty, 1, 1, imgout);
            }
        }
        fclose(imgout);
        return;
    }
    void process(Image c) {
        int ch = c.infoHeader.biHeight, cw = c.infoHeader.biWidth;
        RGB **colrgb = (RGB**)malloc(sizeof(RGB*) * ch);
        for (int i = 0; i < ch; ++i) {
            colrgb[i] = (RGB*)malloc(sizeof(RGB) * cw);
            for (int j = 0; j < cw; ++j) {
                colrgb[i][j].r = (double)c.index[c.pixelIndex[i][j]].r / 255.0;
                colrgb[i][j].g = (double)c.index[c.pixelIndex[i][j]].g / 255.0;
                colrgb[i][j].b = (double)c.index[c.pixelIndex[i][j]].b / 255.0;
            }
        }
        double colr = 0.0, colg = 0.0, colb = 0.0;
        for (int i = 0; i < ch; ++i) {
            for (int j = 0; j < cw; ++j) {
                colr += colrgb[i][j].r;
                colg += colrgb[i][j].g;
                colb += colrgb[i][j].b;
            }
        }
        colr /= ch * cw, colg /= ch * cw, colb /= ch * cw;
        // cerr << colr << " " << colg << " " << colb << " " << ch * cw << endl;
        double sigma_cr = 0.0, sigma_cg = 0.0, sigma_cb = 0.0;
        for (int i = 0; i < ch; ++i) {
            for (int j = 0; j < cw; ++j) {
                sigma_cr += (colrgb[i][j].r - colr) * (colrgb[i][j].r - colr);
                sigma_cg += (colrgb[i][j].g - colg) * (colrgb[i][j].g - colg);
                sigma_cb += (colrgb[i][j].b - colb) * (colrgb[i][j].b - colb);
            }
        }
        sigma_cr /= ch * cw, sigma_cg /= ch * cw, sigma_cb /= ch * cw;
        sigma_cr = sqrt(sigma_cr), sigma_cg = sqrt(sigma_cg), sigma_cb = sqrt(sigma_cb);
        // cerr << sigma_cr << " " << sigma_cg << " " << sigma_cb << endl;
        int h = infoHeader.biHeight, w = infoHeader.biWidth;
        RGB **rgb = (RGB**)malloc(sizeof(RGB*) * h);
        for (int i = 0; i < h; ++i) {
            rgb[i] = (RGB*)malloc(sizeof(RGB) * w);
            for (int j = 0; j < w; ++j) {
                rgb[i][j].r = (double)index[pixelIndex[i][j]].r / 255.0;
                rgb[i][j].g = (double)index[pixelIndex[i][j]].g / 255.0;
                rgb[i][j].b = (double)index[pixelIndex[i][j]].b / 255.0;
            }
        }
        double avgr = 0.0, avgg = 0.0, avgb = 0.0;
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                avgr += rgb[i][j].r;
                avgg += rgb[i][j].g;
                avgb += rgb[i][j].b;
            }
        }
        avgr /= h * w, avgg /= h * w, avgb /= h * w;
        // cerr << avgr << " " << avgg << " " << avgb << " " << h * w << endl;
        double sigma_sr = 0.0, sigma_sg = 0.0, sigma_sb = 0.0;
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                sigma_sr += (rgb[i][j].r - avgr) * (rgb[i][j].r - avgr);
                sigma_sg += (rgb[i][j].g - avgg) * (rgb[i][j].g - avgg);
                sigma_sb += (rgb[i][j].b - avgb) * (rgb[i][j].b - avgb);
            }
        }
        sigma_sr /= h * w, sigma_sg /= h * w, sigma_sb /= h * w;
        sigma_sr = sqrt(sigma_sr), sigma_sg = sqrt(sigma_sg), sigma_sb = sqrt(sigma_sb);
        // cerr << sigma_sr << " " << sigma_sg << " " << sigma_sb << endl;
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                rgb[i][j].r = sigma_cr * (rgb[i][j].r - avgr) / sigma_sr + colr;
                rgb[i][j].g = sigma_cg * (rgb[i][j].g - avgg) / sigma_sg + colg;
                rgb[i][j].b = sigma_cb * (rgb[i][j].b - avgb) / sigma_sb + colb;
                // cerr << i << " " << j << ": " << rgb[i][j].r << " " << rgb[i][j].g << " " << rgb[i][j].b << endl;
            }
        }
        int pos = 0;
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                index[pos].r = (uint8)min(255.0, max(0.0, pow(rgb[i][j].r, 0.7) * 255.0));
                index[pos].g = (uint8)min(255.0, max(0.0, pow(rgb[i][j].g, 0.7) * 255.0));
                index[pos].b = (uint8)min(255.0, max(0.0, pow(rgb[i][j].b, 0.7) * 255.0));
                pixelIndex[i][j] = pos++;
            }
        }
    }
};

int main()
{
    Image color, base;
    base.open("base.bmp");
    color.open("color.bmp");
    base.process(color);
    base.save("result2.bmp");
    return 0;
}