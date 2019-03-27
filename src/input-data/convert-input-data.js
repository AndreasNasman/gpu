exports.execute = () => {
  console.log("\nArcminutes to radians conversion started.");

  const fs = require("fs");

  const arcminutesToRadians = arcminutes =>
    (Number(arcminutes) * Math.PI) / (60 * 180);

  const convertToRadiansJSON = (inputFile, outputName) => {
    const result = [];

    const data = fs.readFileSync(
      `${__dirname}/original/${inputFile}.txt`,
      "utf-8"
    );

    data
      .split("\r\n")
      .splice(1)
      .forEach(row => {
        let [rightAcsension, declination] = row.split(/\s+/);
        declination = arcminutesToRadians(declination);
        rightAcsension = arcminutesToRadians(rightAcsension);

        result.push({
          declination: declination,
          rightAcsension: rightAcsension
        });
      });

    fs.writeFileSync(`${__dirname}/${outputName}.json`, JSON.stringify(result));
    console.log(`'${outputName}' was converted successfully!`);
  };

  convertToRadiansJSON("data_100k_arcmin", "real-galaxies");
  convertToRadiansJSON("flat_100k_arcmin", "synthetic-random-galaxies");

  console.log("Conversion finished.");
};
