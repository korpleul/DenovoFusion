#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

// Function to remove blank lines and trailing spaces from a sequence
void removeBlankLines(std::string& sequence) {
    std::string cleaned;
    for (char ch : sequence) {
        if (!isspace(ch)) {
            cleaned += ch;
        }
    }
    sequence = cleaned;
}

// Function to split sequences and write to output file
void processSequence(const std::string& sequenceName, const std::string& sequence, int maxSequenceLength, std::ofstream& outputFile) {
    int seqLength = sequence.length();
    int currentPosition = 0;
    int sequenceCount = 1;
    std::string subSequence;

    std::string nameCopy = sequenceName.substr(0, sequenceName.find(' ')); // Truncate at first space

    while (currentPosition < seqLength) {
        int remainingLength = seqLength - currentPosition;
        int lengthToCopy = std::min(maxSequenceLength, remainingLength);

        if (lengthToCopy > 0) {
            subSequence = sequence.substr(currentPosition, lengthToCopy);
            removeBlankLines(subSequence);
            if (!subSequence.empty()) {
                outputFile << ">" << nameCopy << "_" << sequenceCount << "\n";
                outputFile << subSequence << "\n";
                sequenceCount++;
            }
        }
        currentPosition += lengthToCopy;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <sequence_length>\n";
        return 1;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];
    int maxSeqLength = std::stoi(argv[3]);

    std::ifstream inputFile(inputFileName);
    std::ofstream outputFile(outputFileName);

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening file\n";
        return 1;
    }

    std::string line;
    std::string sequence;
    std::string sequenceName;

    while (std::getline(inputFile, line)) {
        if (!line.empty() && line.front() == '>') {
            if (!sequenceName.empty()) {
                removeBlankLines(sequence);
                if (!sequence.empty()) {
                    processSequence(sequenceName, sequence, maxSeqLength, outputFile);
                }
                sequence.clear();
            }
            sequenceName = line.substr(1); // Remove '>'
        } else {
            sequence += line;
        }
    }

    if (!sequenceName.empty()) {
        removeBlankLines(sequence);
        if (!sequence.empty()) {
            processSequence(sequenceName, sequence, maxSeqLength, outputFile);
        }
    }

    inputFile.close();
    outputFile.close();

    return 0;
}
