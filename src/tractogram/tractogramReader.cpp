#include "tractogramReader.h"

TractogramReader::TractogramReader() {
	fileName = "";
	fileDescription = "";
	file = NULL;
	threadId = 0;
}

TractogramReader::TractogramReader(std::string _fileName) {
	fileName = "";
	fileDescription = "";
	file = NULL;
	threadId = 0;
	initReader(_fileName);
}

TractogramReader::TractogramReader(const TractogramReader& obj) {
	fileName = obj.fileName;
	fileDescription = obj.fileDescription;
	fileFormat = obj.fileFormat;
	numberOfPoints = obj.numberOfPoints;
	numberOfStreamlines = obj.numberOfStreamlines;
	len = obj.len;
	streamlinePos = obj.streamlinePos;
	file = fopen(fileName.c_str(), "rb+");
}

TractogramReader::~TractogramReader() {
	fileName.erase();
	fileDescription.erase();
	if (file != NULL)
		fclose(file);
	delete[] len;
	delete[] streamlinePos;
	file = NULL;
	len = NULL;
	streamlinePos = NULL;
}



void TractogramReader::copyFrom(const TractogramReader& obj) {
	fileName = obj.fileName;
	fileDescription = obj.fileDescription;
	fileFormat = obj.fileFormat;
	numberOfPoints = obj.numberOfPoints;
	numberOfStreamlines = obj.numberOfStreamlines;
	len = obj.len;
	streamlinePos = obj.streamlinePos;
	file = fopen(fileName.c_str(), "rb+");
}

void TractogramReader::destroyCopy() {
	fileName.erase();
	fileDescription.erase();
	if (file != NULL)
		fclose(file);
	file = NULL;
	len = NULL;
	streamlinePos = NULL;
}


bool TractogramReader::initReader(std::string _fileName) {

	fileName = _fileName;
	file = fopen(fileName.c_str(), "rb+");

	if (file == NULL)
		return false;

	const size_t strLength = 256;
	char dummy[strLength];

	std::string extension = getFileExtension(fileName);

	if (extension == "tck") {

		fileFormat = TCK;
		fileDescription = "mrtrix tracks";

		std::string tmps;
		long pos = 0;
		do {
			fgets(dummy, strLength, file);
			tmps = std::string(dummy);

			size_t column = tmps.find_last_of(":");

			if (column != std::string::npos) {
				std::string key = tmps.substr(0, column);

				if (key == "count")   numberOfStreamlines = std::atoi((tmps.substr(column + 1)).c_str());
				if (key == "file")    pos = std::atoi((tmps.substr(column + 3)).c_str());
			}

		} while (tmps != "END\n");

		if (numberOfStreamlines > 0) {

			numberOfPoints = 0;
			len = new uint32_t[numberOfStreamlines];
			streamlinePos = new long[numberOfStreamlines];
			streamlinePos[0] = pos;
			len[0] = 0;
			fseek(file, pos, SEEK_SET);

			float  tmp[3] = { 0,0,0 };
			size_t ind = 0;

			while (!feof(file)) {
				std::fread(tmp, sizeof(float), 3, file);
				len[ind]++;

				if (std::isnan(tmp[0])) {
					len[ind]--;
					numberOfPoints += len[ind];
					ind++;
					if (ind == numberOfStreamlines)
						break;
					else {
						streamlinePos[ind] = ftell(file);
						len[ind] = 0;
					}
				}
			}

		}

	}

	if (extension == "trk") {
		fileFormat = TRK;
		// TODO: init TRK reader
	}

	if (extension == "vtk") {

		std::string vtkFormat;

		fgets(dummy, strLength, file);                        // vtk version
		fgets(dummy, strLength, file);                        // file description
		fileDescription = std::string(dummy);
		fileDescription = fileDescription.substr(0, fileDescription.size() - 1);
		std::fscanf(file, "%s\n", dummy);                     // ascii or binary
		vtkFormat = std::string(dummy);
		fgets(dummy, strLength, file);                        // always DATASET POLYDATA, we skip to check this for now
		std::fscanf(file, "%*s %zu %*s\n", &numberOfPoints);  // number of points and datatype, we assume datatype is float and skip checking this
		long posData = ftell(file);

		if ((vtkFormat == "ascii") || (vtkFormat == "ASCII")) fileFormat = VTK_ASCII;
		else if ((vtkFormat == "binary") || (vtkFormat == "BINARY")) fileFormat = VTK_BINARY;
		else return false;

		if (numberOfPoints > 0) {

			// Skip points
			if (fileFormat == VTK_BINARY) {
				std::fseek(file, sizeof(float) * numberOfPoints * 3, SEEK_CUR);
				int tmp = std::fgetc(file); if (tmp != '\n') std::ungetc(tmp, file); // Make sure to go end of the line
			}

			if (fileFormat == VTK_ASCII) {
				for (size_t i = 0; i < numberOfPoints; i++)
					fgets(dummy, strLength, file);
			}

			// Get number of streamlines, lengths and streamlinePos
			std::fscanf(file, "LINES %zu %*d\n", &numberOfStreamlines);

			len = new uint32_t[numberOfStreamlines];
			streamlinePos = new long[numberOfStreamlines];
			streamlinePos[0] = posData;

			if (fileFormat == VTK_BINARY) {
				int tmp;
				std::fread(&tmp, sizeof(int), 1, file);
				swapByteOrder(tmp);
				std::fseek(file, tmp * sizeof(int), SEEK_CUR);
				len[0] = tmp;

				for (size_t i = 1; i < numberOfStreamlines; i++) {
					streamlinePos[i] = streamlinePos[i - 1] + long(sizeof(float) * len[i - 1] * 3);
					std::fread(&tmp, sizeof(int), 1, file);
					swapByteOrder(tmp);
					std::fseek(file, tmp * sizeof(int), SEEK_CUR);
					len[i] = tmp;
				}
			}

			if (fileFormat == VTK_ASCII) {
				for (size_t i = 0; i < numberOfStreamlines; i++) {
					std::fscanf(file, "%u ", &len[i]);
					for (uint32_t j = 0; j < (len[i] - 1); j++)
						std::fscanf(file, "%*d ");
					std::fscanf(file, "%*d \n");
				}

				fseek(file, posData, SEEK_SET);
				for (size_t i = 0; i < (numberOfStreamlines - 1); i++) {
					for (size_t j = 0; j < len[i]; j++)
						fgets(dummy, strLength, file);
					streamlinePos[i + 1] = ftell(file);
				}

			}

		}
		else {
			numberOfStreamlines = 0;
			len = NULL;
			streamlinePos = NULL;
		}

	}

	return true;
}

float** TractogramReader::readStreamline(size_t n) {

	float** points;
	points = new float* [len[n]];


	fseek(file, streamlinePos[n], SEEK_SET);

	if (fileFormat == TCK) {
		float tmp;
		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, file);
				points[i][j] = tmp;
			}
		}
	}

	if (fileFormat == TRK) {
		// TODO 
	}

	if (fileFormat == VTK_BINARY) {
		float tmp;
		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, file);
				swapByteOrder(tmp);
				points[i][j] = tmp;
			}
		}
	}

	if (fileFormat == VTK_ASCII) {

		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];
			std::fscanf(file, "%f %f %f\n", &points[i][0], &points[i][1], &points[i][2]);
		}

	}

	return points;

}

void TractogramReader::readPoint(size_t n, uint32_t l, float* point) {

	if (fileFormat == TCK) {
		fseek(file, streamlinePos[n] + sizeof(float) * 3 * l, SEEK_SET);
		float tmp;
		for (int j = 0; j < 3; j++) {
			std::fread(&tmp, sizeof(float), 1, file);
			point[j] = tmp;
		}
	}

	if (fileFormat == TRK) {
		// TODO 
	}

	if (fileFormat == VTK_BINARY) {
		fseek(file, streamlinePos[n] + sizeof(float) * 3 * l, SEEK_SET);
		float tmp;
		for (int j = 0; j < 3; j++) {
			std::fread(&tmp, sizeof(float), 1, file);
			swapByteOrder(tmp);
			point[j] = tmp;
		}
	}

	// For ASCII files, reading points is not efficient
	if (fileFormat == VTK_ASCII) {
		fseek(file, streamlinePos[n], SEEK_SET);
		const size_t strLength = 256;
		char dummy[strLength];

		for (uint32_t i = 0; i < (l - 1); i++) {
			fgets(dummy, strLength, file);
		}
		std::fscanf(file, "%f %f %f\n", &point[0], &point[1], &point[2]);
	}

}


std::vector<Point> TractogramReader::readVectorizedStreamline(size_t n) {

	std::vector<Point> points;
	points.reserve(len[n]);

	fseek(file, streamlinePos[n], SEEK_SET);

	if (fileFormat == TCK) {
		for (uint32_t i = 0; i < len[n]; i++) {
			Point p;
			std::fread(&p.x, sizeof(float), 1, file);
			std::fread(&p.y, sizeof(float), 1, file);
			std::fread(&p.z, sizeof(float), 1, file);
			points.push_back(p);
		}
	}


	if (fileFormat == TRK) {
		// TODO 
	}

	if (fileFormat == VTK_BINARY) {
		for (uint32_t i = 0; i < len[n]; i++) {
			Point p;
			std::fread(&p.x, sizeof(float), 1, file); swapByteOrder(p.x);
			std::fread(&p.y, sizeof(float), 1, file); swapByteOrder(p.y);
			std::fread(&p.z, sizeof(float), 1, file); swapByteOrder(p.z);
			points.push_back(p);
		}
	}

	if (fileFormat == VTK_ASCII) {

		for (uint32_t i = 0; i < len[n]; i++) {
			Point p;
			std::fscanf(file, "%f %f %f\n", &p.x, &p.y, &p.z);
			swapByteOrder(p.x);
			swapByteOrder(p.y);
			swapByteOrder(p.z);
			points.push_back(p);
		}

	}

	return points;

}
