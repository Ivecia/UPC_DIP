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

inline int** initialization(int k)
{
    int **mod = (int**)malloc(sizeof(int*) * k);
    for (int i = 0; i < k; ++i) {
        mod[i] = (int*)malloc(sizeof(int) * k);
        for (int j = 0; j < k; ++j) mod[i][j] = 0;
    }
    return mod;
}

inline uint8** normalization(int** ori, int h, int w)
{
    int l = 0x3f3f3f3f, r = -0x3f3f3f3f;
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            l = min(l, ori[i][j]);
            r = max(r, ori[i][j]);
        }
    }
    uint8 **ret = (uint8**)malloc(sizeof(uint8*) * h);
    for (int i = 0; i < h; ++i) {
        ret[i] = (uint8*)malloc(sizeof(uint8) * w);
        for (int j = 0; j < w; ++j) {
            if (l == r) ret[i][j] = (uint8)max(0, min(ori[i][j], 255));
            else ret[i][j] = (uint8)((double)(ori[i][j] - l) / (r - l) * 255);
        }
    }
    free(ori);
    return ret;
}

inline int** sobelOperatorHorizontal()
{
    int **mod = initialization(3);
    mod[0][0] = mod[2][0] = -1;
    mod[0][2] = mod[2][2] = 1;
    mod[1][0] = -2, mod[1][2] = 2;
    return mod;
}

inline int** sobelOperatorVertical()
{
    int **mod = initialization(3);
    mod[0][0] = mod[0][2] = -1;
    mod[2][0] = mod[2][2] = 1;
    mod[0][1] = -2, mod[2][1] = 2;
    return mod;
}

inline int** laplacianOperator4()
{
    int **mod = initialization(3);
    mod[0][1] = mod[1][0] = mod[1][2] = mod[2][1] = -1;
    mod[2][2] = 4;
    return mod;
}

inline int** laplacianOperator8()
{
    int **mod = initialization(3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            mod[i][j] = -1;
        }
    }
    mod[2][2] = 8;
    return mod;
}

uint8** addition(uint8** l, uint8 **r, int height, int width)
{
    uint8 **ret = (uint8**)malloc(sizeof(uint8*) * height);
    for (int i = 0; i < height; ++i) {
        ret[i] = (uint8*)malloc(sizeof(uint8) * width);
        for (int j = 0; j < width; ++j) {
            ret[i][j] = max(0, min((int)l[i][j] + r[i][j], 255));
        }
    }
    return ret;
}

uint8** multiplication(uint8** l, uint8 **r, int height, int width)
{
    int **tmp = (int**)malloc(sizeof(int*) * height);
    for (int i = 0; i < height; ++i) {
        tmp[i] = (int*)malloc(sizeof(int) * width);
        for (int j = 0; j < width; ++j) {
            tmp[i][j] = (int)l[i][j] * r[i][j];
        }
    }
    return normalization(tmp, height, width);
}

uint8** mix(uint8** l, int pl, uint8 **r, int pr, int height, int width)
{
    int p = pl + pr;
    uint8 **ret = (uint8**)malloc(sizeof(uint8*) * height);
    for (int i = 0; i < height; ++i) {
        ret[i] = (uint8*)malloc(sizeof(uint8) * width);
        for (int j = 0; j < width; ++j) {
            ret[i][j] = (pl * l[i][j] + pr * r[i][j] + (p >> 1)) / p;
        }
    }
    return ret;
}

uint8** gamma(uint8** ori, double gamma, int height, int width)
{
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double tmp = (double)ori[i][j] / 255.0;
            tmp = pow(tmp, gamma);
            ori[i][j] = (uint8)(255.0 * tmp + 0.5);
        }
    }
    return ori;
}

uint8** medianFiltering(uint8** pixelIndex, int k, int height, int width)
{
    uint8 *tmp = (uint8*)malloc(sizeof(uint8) * k);
    uint8 **ret = (uint8**)malloc(sizeof(uint8*) * height);
    for (int i = 0; i < height; ++i) {
        ret[i] = (uint8*)malloc(sizeof(uint8) * width);
        for (int j = 0; j < width; ++j) {
            int num = 0;
            for (int c = max(0, i - (k >> 1)); c < min(height, i + (k >> 1) + 1); ++c) {
                for (int d = max(0, j - (k >> 1)); d < min(width, j + (k >> 1) + 1); ++d) {
                    tmp[num++] = pixelIndex[c][d];
                }
            }
            sort(tmp, tmp + num);
            ret[i][j] = tmp[(num >> 1) + 1];
        }
    }
    free(tmp);
    return ret;
}

uint8** nonlinearSharpening(uint8** pixelIndex, int k, int **mod, int height, int width)
{
    uint8 **ret = (uint8**)malloc(sizeof(uint8*) * height);
    for (int i = 0; i < height; ++i) {
        ret[i] = (uint8*)malloc(sizeof(uint8) * width);
        for (int j = 0; j < width; ++j) {
            int val = 0;
            for (int c = max(0, i - (k >> 1)); c < min(height, i + (k >> 1) + 1); ++c) {
                for (int d = max(0, j - (k >> 1)); d < min(width, j + (k >> 1) + 1); ++d) {
                    val += mod[c - (i - (k >> 1))][d - (j - (k >> 1))] * pixelIndex[c][d];
                }
            }
            ret[i][j] = (uint8)max(0, min(val + pixelIndex[i][j], 255));
        }
    }
    return ret;
}

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
        fread(&alpha, 1, 1, imgin);
    }
    inline void write(FILE *imgout) {
        fwrite(&b, 1, 1, imgout);
        fwrite(&g, 1, 1, imgout);
        fwrite(&r, 1, 1, imgout);
        fwrite(&alpha, 1, 1, imgout);
    }
};

class Image {
private:
    BitMapFileHeader fileHeader;
    BitMapInfoHeader infoHeader;
    BitMapColorIndex *index;
    uint8 **pixelIndex;
public:
    void open() {
        FILE *imgin = fopen("test.bmp", "rb");
        if (imgin == 0) errorMessage(ImageNotExists);
        fileHeader.read(imgin);
        infoHeader.read(imgin);
        index = (BitMapColorIndex*)malloc(sizeof(BitMapColorIndex) * infoHeader.biClrUsed);
        for (int i = 0; i < infoHeader.biClrUsed; ++i) {
            index[i].read(imgin);
        }
        pixelIndex = (uint8**)malloc(sizeof(uint8*) * infoHeader.biHeight);
        if (infoHeader.biCompression == 0) {
            uint8 num;
            for (int i = infoHeader.biHeight - 1; i >= 0; --i) {
                pixelIndex[i] = (uint8*)malloc(sizeof(uint8) * infoHeader.biWidth);
                int lineByte = (infoHeader.biWidth * infoHeader.biBitCount / 8 + 3) / 4 * 4;
                for (int j = 0; j < lineByte; ++j) {
                    fread(&num, 1, 1, imgin);
                    if (j < infoHeader.biWidth) {
                        pixelIndex[i][j] = num;
                    }
                }
            }
        } else if (infoHeader.biCompression == 1) {
            uint8 opt, num;
            for (int i = infoHeader.biHeight - 1; i >= 0; --i) {
                pixelIndex[i] = (uint8*)malloc(sizeof(uint8) * infoHeader.biWidth);
                int j = 0;
                while (1) {
                    fread(&opt, 1, 1, imgin);
                    fread(&num, 1, 1, imgin);
                    if (opt == 0) {
                        if (num == 0) break;
                        else errorMessage(UnknownError);
                    } else {
                        for (int k = 0; k < opt; ++k, ++j) {
                            pixelIndex[i][j] = num;
                        }
                    }
                }
            }
            fread(&opt, 1, 1, imgin);
            fread(&num, 1, 1, imgin);
            if (opt != 0 || num != 1) errorMessage(BinaryFileError);
        } else errorMessage(UnrecognizedCompressionType);
        fclose(imgin);
        return;
    }
    void save() {
        FILE *imgout = fopen("result.bmp", "wb");
        infoHeader.biCompression = 0;
        infoHeader.biBitCount = 8;
        infoHeader.biSizeImage = 0;
        infoHeader.biXPelsPerMeter = 0;
        infoHeader.biYPelsPerMeter = 0;
        fileHeader.write(imgout);
        infoHeader.write(imgout);
        for (int i = 0; i < infoHeader.biClrUsed; ++i) {
            index[i].write(imgout);
        }
        uint8 empty = 0;
        for (int i = infoHeader.biHeight - 1; i >= 0; --i) {
            int lineByte = (infoHeader.biWidth * infoHeader.biBitCount / 8 + 3) / 4 * 4;
            for (int j = 0; j < infoHeader.biWidth; ++j) {
                fwrite(&pixelIndex[i][j], 1, 1, imgout);
            }
            for (int j = infoHeader.biWidth; j < lineByte; ++j) {
                fwrite(&empty, 1, 1, imgout);
            }
        }
        fclose(imgout);
        return;
    }
    void process() {
        int h = infoHeader.biHeight, w = infoHeader.biWidth;
        /*
        uint8 **c = nonlinearSharpening(pixelIndex, 3, laplacianOperator8(), h, w);
        uint8 **dx = nonlinearSharpening(pixelIndex, 3, sobelOperatorHorizontal(), h, w);
        uint8 **dy = nonlinearSharpening(pixelIndex, 3, sobelOperatorVertical(), h, w);
        uint8 **d = mix(dx, 1, dy, 1, h, w);
        uint8 **e = medianFiltering(d, 3, h, w);
        uint8 **f = multiplication(c, e, h, w);
        uint8 **g = addition(pixelIndex, f, h, w);
        pixelIndex = gamma(g, 0.5, h, w);
        */
        uint8 **c = gamma(pixelIndex, 0.85, h, w);
        uint8 **dx = nonlinearSharpening(c, 3, sobelOperatorHorizontal(), h, w);
        uint8 **dy = nonlinearSharpening(c, 3, sobelOperatorVertical(), h, w);
        uint8 **d = mix(dx, 1, dy, 1, h, w);
        uint8 **e = nonlinearSharpening(pixelIndex, 3, laplacianOperator8(), h, w);
        uint8 **f = multiplication(d, e, h, w);
        uint8 **m = addition(f, pixelIndex, h, w);
        pixelIndex = gamma(m, 0.85, h, w);
    }
};

int main()
{
    Image img;
    img.open();
    img.process();
    img.save();
    return 0;
}