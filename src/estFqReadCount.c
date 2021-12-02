/*
GZIP file format specification version 4.3
https://www.ietf.org/rfc/rfc1952.txt

Code adapted from Heng Li's samtools
Written by Chee Hong Wong @cheehongsg
*/
#include <zlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static int G_Verbose = 0;
double G_t_real;

/*********
 * Timer *
 *********/

double cputime()
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

/*************
 * Estimator *
 *************/

off_t getFileSize(const char *szFile)
{
    int fd;
    struct stat statbuf;
    off_t filesize = 0;

    fd = open(szFile, O_RDONLY, S_IRUSR | S_IRGRP);
    if (fd == -1)
    {
        fprintf(stderr, "[E::%s] fail to open '%s' for reading: %s\n", __func__, szFile, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (fstat(fd, &statbuf) == -1)
    {
        fprintf(stderr, "[E::%s] fail to stat '%s': %s\n", __func__, szFile, strerror(errno));
        exit(EXIT_FAILURE);
    }

    filesize = statbuf.st_size;

    if (close(fd) == -1)
    {
        fprintf(stderr, "[E::%s] fail to close '%s': %s\n", __func__, szFile, strerror(errno));
        exit(EXIT_FAILURE);
    }

    return filesize;
}

void estimateReadCount(const char *szFile, int subsamplePercentage, long reportBlock)
{
    gzFile fp;
    kseq_t *seq;
    int l, i;
    off_t fileSize = getFileSize(szFile);
    z_off_t compressedSize = 0;
    z_off_t prevCompressedSize = 0;
    long long uncompressedSize = 0;
    long long numFastqRecords = 0;
    float compressionRatio = 1.0;
    long long averageFastqRecordSize = 1;
    double estimatedTotalRecords = 0.0;

    long long compressedSizeThreshold = (long long)(((double) fileSize / 100.0f) * subsamplePercentage);
    compressedSizeThreshold += 16384;

    fp = gzopen(szFile, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {
        if (!seq->is_fastq)
        {
            break;
        }

        prevCompressedSize = compressedSize;
        compressedSize = gzoffset(fp);
        if (compressedSize>compressedSizeThreshold) {
            compressedSize = prevCompressedSize;
            break;
        }

        uncompressedSize += (seq->name.l + 1 + l + 1 + 2 + l + 1);
        if (seq->comment.l>0) {
            uncompressedSize += (1 + seq->comment.l);
        }
        numFastqRecords++;

        if (G_Verbose>=3)
        {
            if (0 == (numFastqRecords % reportBlock))
            {
                // est. number of record = <fileSize> / <compressed size> * <uncompressed size> / <average fastq record size>
                
                compressionRatio = (double)compressedSize / (double)uncompressedSize;
                averageFastqRecordSize = uncompressedSize / numFastqRecords;
                estimatedTotalRecords = (double)fileSize / (double)averageFastqRecordSize / (double)compressionRatio;

                fprintf(stderr, "[E::%s] file size: %llu compressed size: %ld (%0.4f) uncompressed size: %llu fastqCount: %llu compression ratio: %.4f record size: %llu est. records count: %.0f\n",
                        __func__,
                        fileSize, compressedSize, (double)compressedSize / (double)fileSize, uncompressedSize, numFastqRecords,
                        compressionRatio, averageFastqRecordSize, estimatedTotalRecords);
            }
        }
    }
    if (l >= 0)
    {
        if (!seq->is_fastq)
        {
            fprintf(stderr, "[E::%s] '%s' is not a fastq file.\n", __func__, szFile);
            kseq_destroy(seq);
            gzclose(fp);
            exit(EXIT_FAILURE);
        }
    } else {
        if (-2==l) {
            fprintf(stderr, "[E::%s] Truncated fastq record '%s' in '%s'\n", __func__, seq->name.s, szFile);
        } else if (-3==l) {
            fprintf(stderr, "[E::%s] Stream error in '%s', last successful record #'%llu'\n", __func__, szFile, numFastqRecords);
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    if (l < -1)
    {
        exit(EXIT_FAILURE);
    }
    else
    {
        if (0 == numFastqRecords)
        {
            // TODO: inform user
            return;
        }
        compressionRatio = (double)compressedSize / (double)uncompressedSize;
        averageFastqRecordSize = uncompressedSize / numFastqRecords;
        estimatedTotalRecords = (double)fileSize / (double)averageFastqRecordSize / (double)compressionRatio;
        if (G_Verbose>=3)
        {
            fprintf(stderr, "[E::%s] file size: %llu compressed size: %ld uncompressed size: %llu fastqCount: %llu compression ratio: %.4f record size: %llu est. records count: %.0f\n",
                    __func__,
                    fileSize, compressedSize, uncompressedSize, numFastqRecords,
                    compressionRatio, averageFastqRecordSize, estimatedTotalRecords);
        } 
        else if (G_Verbose>=1)
        {
            fprintf(stderr, "[M::%s] file: %s file size: %llu est. records count: %.0f\n", __func__, szFile, fileSize, estimatedTotalRecords);
        }
        else
        {
            fprintf(stdout, "%.0f\n", estimatedTotalRecords);
        }
    }
}

int main(int argc, char *argv[])
{
    gzFile fp;
    kseq_t *seq;
    int l, i;
    long long llTotalBases = 0;
    int c;
    int subsamplePercentage = 1;
    long blockSize = 100000;
    double t_diff;

    G_t_real = realtime();

    while ((c = getopt(argc, argv, "v:s:b:")) >= 0)
    {
        switch (c)
        {
        case 'v':
            G_Verbose = atoi(optarg);
            G_Verbose = G_Verbose>0 ? G_Verbose : 0;
            break;
        case 's':
            subsamplePercentage = atoi(optarg);
            subsamplePercentage = subsamplePercentage>0 ? subsamplePercentage : 1;
            break;
        case 'b':
            blockSize = atol(optarg);
            blockSize = blockSize>0 ? blockSize : 100000;
            break;
        default:
            return 1;
        }
    }

    if (optind + 1 > argc)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: %s [options] <fastqFile>\n", argv[0]);
		fprintf(stderr, "       -s INT     percentage to subsample to estimate reads count (1 to 100) [%d]\n", 1);
        fprintf(stderr, "\n");
        return EXIT_FAILURE;
    }

    estimateReadCount(argv[optind], subsamplePercentage, blockSize);

    if (G_Verbose>=1)
    {
        fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - G_t_real, cputime());
    }
    return 0;
}
