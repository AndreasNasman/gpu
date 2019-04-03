exports.execute = () => {
  console.log("\x1b[33m%s\x1b[0m", "Building of histograms started.");

  const fs = require("fs");

  const realGalaxies = require("../input-data/real-galaxies.json");
  const randomGalaxies = require("../input-data/random-galaxies.json");

  const COVERAGE = 180; // degrees
  const BIN_WIDTH = 0.25; // degress
  const NUMBER_OF_BINS = COVERAGE * (1 / BIN_WIDTH);

  const buildHistogram = (firstGalaxySet, secondGalaxySet, outputName) => {
    const bins = Array(NUMBER_OF_BINS).fill(0);

    const angleBetweenTwoGalaxies = (firstGalaxy, secondGalaxy) =>
      Math.acos(
        Math.sin(firstGalaxy.declination) * Math.sin(secondGalaxy.declination) +
          Math.cos(firstGalaxy.declination) *
            Math.cos(secondGalaxy.declination) *
            Math.cos(firstGalaxy.rightAscension - secondGalaxy.rightAscension)
      );

    const radiansToDegrees = radians => radians * (180 / Math.PI);

    const updateBins = angle => {
      if (isNaN(angle)) return;
      const index = Math.floor(radiansToDegrees(angle) / BIN_WIDTH);
      bins[index] += 1;
    };

    const SIZE = 10000;
    if (!secondGalaxySet) {
      firstGalaxySet.slice(0, SIZE).forEach((firstGalaxy, index) => {
        firstGalaxySet.slice(index, SIZE).forEach(secondGalaxy => {
          const angle = angleBetweenTwoGalaxies(firstGalaxy, secondGalaxy);
          updateBins(angle);
        });
      });
    } else {
      firstGalaxySet.slice(0, SIZE).forEach(firstGalaxy => {
        secondGalaxySet.slice(0, SIZE).forEach(secondGalaxy => {
          const angle = angleBetweenTwoGalaxies(firstGalaxy, secondGalaxy);
          updateBins(angle);
        });
      });
    }

    fs.writeFileSync(`${__dirname}/${outputName}.json`, JSON.stringify(bins));
    console.log(`Histogram '${outputName}' built!`);
  };

  buildHistogram(realGalaxies, null, "DD");
  buildHistogram(realGalaxies, randomGalaxies, "DR");
  buildHistogram(randomGalaxies, null, "RR");

  console.log("\x1b[33m%s\x1b[0m", "Building of histograms started.\n");
};
