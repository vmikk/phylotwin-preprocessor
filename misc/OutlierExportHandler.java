import de.lmu.ifi.dbs.elki.result.Result;
import de.lmu.ifi.dbs.elki.result.ResultHandler;
import de.lmu.ifi.dbs.elki.result.ResultUtil;
import de.lmu.ifi.dbs.elki.result.outlier.OutlierResult;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.ids.DBIDIter;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.model.Model;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.AbstractParameterizer;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.OptionID;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.Parameterization;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameters.FileParameter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class OutlierExportHandler implements ResultHandler {
    private File outputFile;
    
    public OutlierExportHandler(File output) {
        this.outputFile = output;
    }
    
    @Override
    public void processNewResult(Result newResult) {
        Database db = ResultUtil.findDatabase(newResult);
        ArrayList<Clustering<?>> clusterings = ResultUtil.filterResults(newResult, Clustering.class);
        
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            for (Clustering<?> c : clusterings) {
                // Find noise cluster (outliers)
                for (Cluster<?> cluster : c.getAllClusters()) {
                    if (cluster.isNoise()) {
                        writer.write("# Outliers:\n");
                        for (DBIDIter iter = cluster.getIDs().iter(); iter.valid(); iter.advance()) {
                            // Write the object's data
                            writer.write(db.getObject(iter).toString());
                            writer.newLine();
                        }
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Error writing outliers to file: " + e.getMessage(), e);
        }
    }
    
    public static class Parameterizer extends AbstractParameterizer {
        public static final OptionID OUTPUT_ID = new OptionID("outlierexport.output",
                "File to write the outliers to.");
        
        private File outputFile;
        
        @Override
        protected void makeOptions(Parameterization config) {
            super.makeOptions(config);
            FileParameter param = new FileParameter(OUTPUT_ID, FileParameter.FileType.OUTPUT_FILE);
            if (config.grab(param)) {
                outputFile = param.getValue();
            }
        }
        
        @Override
        protected OutlierExportHandler makeInstance() {
            return new OutlierExportHandler(outputFile);
        }
    }
}