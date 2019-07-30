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

    const updateBins = (angle, amount = 1) => {
      const index = Math.floor(radiansToDegrees(angle) / BIN_WIDTH);
      bins[index] += amount;
    };

    firstGalaxySet.forEach((firstGalaxy, index) => {
      if (secondGalaxySet) {
        secondGalaxySet.forEach(secondGalaxy => {
          const angle = angleBetweenTwoGalaxies(firstGalaxy, secondGalaxy);
          updateBins(angle);
        });
      } else {
        const angle = 0;
        updateBins(angle);

        firstGalaxySet.slice(index + 1).forEach(secondGalaxy => {
          const angle = angleBetweenTwoGalaxies(firstGalaxy, secondGalaxy);
          updateBins(angle, 2);
        });
      }
    });

    fs.writeFileSync(`${__dirname}/${outputName}.json`, JSON.stringify(bins));
    console.log(`Histogram '${outputName}' built!`);
  };

  buildHistogram(realGalaxies, null, "DD");
  buildHistogram(realGalaxies, randomGalaxies, "DR");
  buildHistogram(randomGalaxies, null, "RR");

  console.log("\x1b[33m%s\x1b[0m", "Building of histograms finished.\n");
};
