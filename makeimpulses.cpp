#include <math.h>
#include <stdint.h>

#include <fstream>
#include <sstream>

#include "HrtfData.h"

#define SIG(a,b,c,d) ((a) | ((b) << 8) | ((c) << 16) | ((d) << 24))

typedef struct tWAVEFORMATEX {
  uint16_t wFormatTag;
  uint16_t nChannels;
  uint32_t nSamplesPerSec;
  uint32_t nAvgBytesPerSec;
  uint16_t nBlockAlign;
  uint16_t wBitsPerSample;
  uint16_t cbSize;
} WAVEFORMATEX, *PWAVEFORMATEX, *NPWAVEFORMATEX, *LPWAVEFORMATEX;

#define WAVE_FORMAT_FLOAT_NORM (3)

typedef struct speakerPosition {
    const char *name;
	float elevation;
	float azimuth;
	float distance;
} speakerPosition;

#define DEGREES(x) ((x)*M_PI / 180.0)

static const speakerPosition speakerPositions[12] = {
	{ .name = "FL", .elevation = DEGREES(0.0), .azimuth = DEGREES(-30.0), .distance = 1.0 },
	{ .name = "FR", .elevation = DEGREES(0.0), .azimuth = DEGREES(+30.0), .distance = 1.0 },
	{ .name = "FC", .elevation = DEGREES(0.0), .azimuth = DEGREES(0.0), .distance = 1.0 },
	{ .name = "LFE", .elevation = DEGREES(0.0), .azimuth = DEGREES(0.0), .distance = 1.0 },
	{ .name = "BL", .elevation = DEGREES(0.0), .azimuth = DEGREES(-135.0), .distance = 1.0 },
	{ .name = "BR", .elevation = DEGREES(0.0), .azimuth = DEGREES(+135.0), .distance = 1.0 },
	{ .name = "SL", .elevation = DEGREES(0.0), .azimuth = DEGREES(-90.0), .distance = 1.0 },
	{ .name = "SR", .elevation = DEGREES(0.0), .azimuth = DEGREES(+90.0), .distance = 1.0 },
	{ .name = "TFL", .elevation = DEGREES(+45.0), .azimuth = DEGREES(-30.0), .distance = 1.0 },
	{ .name = "TFR", .elevation = DEGREES(+45.0), .azimuth = DEGREES(+30.0), .distance = 1.0 },
	{ .name = "TBL", .elevation = DEGREES(+45.0), .azimuth = DEGREES(-135.0), .distance = 1.0 },
	{ .name = "TBR", .elevation = DEGREES(+45.0), .azimuth = DEGREES(+135.0), .distance = 1.0 }
};

void cblas_scopy(size_t count, const float *src, size_t src_stride, float *dst, size_t dst_stride) {
    for (size_t i = 0; i < count; ++i) {
        *dst = *src;
        src += src_stride;
        dst += dst_stride;
    }
}

int main(void) {
    try {
        std::ifstream file("SADIE_D02-96000.mhr", std::fstream::binary);

        if(!file.is_open()) {
            throw std::logic_error("Cannot open file.");
        }

        HrtfData data(file);

        file.close();

        const float *impulseData = NULL;
        double sampleRateOfSource = 0;
        int sampleCount = 0;

        sampleRateOfSource = data.get_sample_rate();

        uint32_t sampleCountExact = data.get_response_length();

        sampleCount = sampleCountExact + ((data.get_longest_delay() + 2) >> 2);

        for(size_t i = 0; i < 12; ++i) {
            std::vector<float> hrtfData(sampleCount * 2, 0.0);

            const speakerPosition &speaker = speakerPositions[i];
            DirectionData hrtfLeft;
            DirectionData hrtfRight;

            data.get_direction_data(speaker.elevation, speaker.azimuth, speaker.distance, hrtfLeft, hrtfRight);

            cblas_scopy(sampleCountExact, &hrtfLeft.impulse_response[0], 1, &hrtfData[((hrtfLeft.delay + 2) >> 2) * 2], 2);
            cblas_scopy(sampleCountExact, &hrtfRight.impulse_response[0], 1, &hrtfData[((hrtfLeft.delay + 2) >> 2) * 2 + 1], 2);

            uint32_t srate = (uint32_t)(int)(sampleRateOfSource + 0.5);

            WAVEFORMATEX wfx = {
                .wFormatTag      = WAVE_FORMAT_FLOAT_NORM,
                .nChannels       = 2,
                .nSamplesPerSec  = srate,
                .nAvgBytesPerSec = srate * 2 * 4,
                .nBlockAlign     = 2 * 4,
                .wBitsPerSample  = 32,
                .cbSize          = 0
            };

            uint32_t signature = SIG('R', 'I', 'F', 'F');
            uint32_t wavesig = SIG('W', 'A', 'V', 'E');
            uint32_t fmtsig = SIG('f', 'm', 't', ' ');
            uint32_t datasig = SIG('d', 'a', 't', 'a');
            uint32_t fmtsize = sizeof(wfx);
            uint32_t datasize = sampleCount * 2 * 4;
            uint32_t riffsize = sizeof(signature) + sizeof(riffsize) + sizeof(wavesig) + sizeof(fmtsig) + sizeof(fmtsize) + fmtsize + sizeof(datasig) + sizeof(datasize) + datasize;

            std::ostringstream oss;
            oss << "impulse_" << speaker.name << "_" << srate << ".wav";
            std::ofstream ofile(oss.str(), std::fstream::binary);

            if(!ofile.is_open()) {
                throw std::logic_error("Unable to open output file.");
            }

            ofile.write((const char *)&signature, sizeof(signature));
            ofile.write((const char *)&riffsize, sizeof(riffsize));
            ofile.write((const char *)&wavesig, sizeof(wavesig));
            ofile.write((const char *)&fmtsig, sizeof(fmtsig));
            ofile.write((const char *)&fmtsize, sizeof(fmtsize));
            ofile.write((const char *)&wfx, fmtsize);
            ofile.write((const char *)&datasig, sizeof(datasig));
            ofile.write((const char *)&datasize, sizeof(datasize));
            ofile.write((const char *)&hrtfData[0], datasize);

            ofile.close();
       }
    } catch(std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}