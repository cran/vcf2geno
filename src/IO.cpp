#include "IO.h"

//static method
AbstractFileReader* AbstractFileReader::open(const char* fileName){
    AbstractFileReader * fr = NULL;
    if (!fileName || fileName[0] == '\0') {
        REprintf("Empty file name.\n");
        return fr;
    }

    switch(AbstractFileReader::checkFileType(fileName)) {
    case PLAIN:
        fr = new PlainFileReader(fileName);
        break;
    case GZIP:
        fr = new GzipFileReader(fileName);
        break;
#ifdef HAS_BZIP2
      case BZIP2:
        fr = new Bzip2FileReader(fileName);
        break;
#endif
      default:
        REprintf("Cannot detect file type (does it exist?!)\n");
        break;
    }
    return fr;
}
// static method
void AbstractFileReader::close(AbstractFileReader** f) {
    ////assert(f && *f);
    (*f)->close();
    delete (*f);
    *f = NULL;
}

// check header for known file type
AbstractFileReader::FileType AbstractFileReader::checkFileType(const char* fileName){
    // treat stdin as plain text file
    if (strncmp(fileName, "-", 1) == 0) {
        return PLAIN;
    }
    // read file header    
    FILE* fp = fopen(fileName, "rb");
    if (!fp) return UNKNOWN;
    unsigned char header[2]={0,0};
    int n = fread(header, sizeof(char), 2, fp);
    fclose(fp);
    // check file types
    // magic numbers
    const int gz_magic[2] = {0x1f, 0x8b}; /* gzip magic header */
    if ( n >= 2 && header[0] == gz_magic[0] && header[1] == gz_magic[1]) {
        return GZIP;
    }
#ifdef HAS_BZIP2
    const int bzip2_magic[2] = {'B', 'Z'}; /* bzip2 magic header */
    if ( n >= 2 && header[0] == bzip2_magic[0] && header[1] == bzip2_magic[1]) {
        return BZIP2;
    }
#endif
    return PLAIN;
    /* // check the characters fall into visible ASCII range */
    /* if ( header[0] >= 0x20 /\* space *\/ && */
    /*      header[0] <  0x7f /\* DEL *\/   && */
    /*      header[1] >= 0x20 /\* space *\/ && */
    /*      header[1] <  0x7f /\* DEL *\/) { */
    /*     return PLAIN; */
    /* }  */
    /* return UNKNOWN; */
}

AbstractFileWriter::~AbstractFileWriter() {
#ifdef IO_DEBUG
    REprintf("AbstractFileWriter desc()\n");
#endif
}


int BGZipFileWriter::open(const char* fn, bool append){
    if (append) 
        REprintf("Gzip does not support appending.\n");
    this->fp = bgzf_open(fn, "w");
    if (!this->fp) {
        REprintf("ERROR: Cannot open %s for write\n", fn);
        return -1;
    }
    return 0;
}
void BGZipFileWriter::close(){
    if (this->fp) {
        bgzf_close(this->fp);
        this->fp = NULL;
    }
}
int BGZipFileWriter::write(const char* s) {
        return bgzf_write(this->fp, s, strlen(s));
}
int BGZipFileWriter::writeLine(const char* s) {
    int ret = bgzf_write(this->fp, s, strlen(s));
    ret += bgzf_write(this->fp, "\n", 1);
        return (ret);
}

