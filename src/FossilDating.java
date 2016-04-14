import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by alexei on 31/07/15.
 */
public class FossilDating {

    public static void main(String[] args) throws IOException {

        String filename = args[0];

        double upper = 160.0;
        if (args.length > 1) {
            upper = Double.parseDouble(args[1]);
        }
        List<String> taxa = new ArrayList<>();
        if (args.length > 2) {
             taxa = readTaxa(args[2]);
        }

        String idString = null;
        if (args.length > 3) {
            idString = args[3];
        }

        System.out.println("Processing file " + filename);

        String fileStem = filename.substring(0,filename.indexOf("."));

        System.out.println("File stem for output is " + fileStem);

        List<String> lines = new ArrayList<>();


        BufferedReader reader = new BufferedReader(new FileReader(filename));
        Map<String,String> sdMap = new HashMap<>();

        String line = reader.readLine();
        while (line != null) {

            System.out.println(line);

            if (line.matches(".*SamplingDate.*")) {
                String taxon = taxon(line);
                System.out.println("Found sampling dates for taxon " + taxon);
                sdMap.put(taxon, line);
            }

            lines.add(line);
            line = reader.readLine();
        }
        reader.close();

        if (taxa.size() == 0) {
            taxa.addAll(sdMap.keySet());
            if (idString == null) idString = "idref";
        }
        if (idString == null) idString = "id";

        for (String taxon : taxa) {
            String newFile = fileStem + "_" + taxon + ".xml";
            PrintWriter writer = new PrintWriter(new FileWriter(newFile));
            for (String l : lines) {
                if (l.matches(".*SamplingDate.*")) {
                    if (taxon(l).equals(taxon)) {
                        l = newLowerUpper(l,0,upper);
                    }
                }
                writer.write(l);
                writer.write("\n");
                if (l.matches(".*INSERT SINGLE SAMPLINGDATE.*")) {
                    writer.write("    <operator spec='SampledNodeDateRandomWalkerForZeroBranchSATrees' windowSize=\"2.\"  tree=\"@tree\" weight=\"0.5\" useWindowSizeWithSamplingDates=\"true\">\n" +
                            "        <taxonset spec=\"TaxonSet\">\n" +
                            "            <taxon " + idString + "=\"" + taxon + "\" spec=\"beast.evolution.alignment.Taxon\"/>\n" +
                            "        </taxonset>\n" +
                            "        <samplingDates id=\"samplingDate\" spec=\"beast.evolution.tree.SamplingDate\" taxon=\"@" + taxon + "\" lower=\"0.0\" upper=\"" + upper + "\"/>\n" +
                            "    </operator>\n");
                }
            }
            writer.flush();
            writer.close();
        }
    }

    private static List<String> readTaxa(String filename) throws IOException {

        System.out.println("Reading taxa from file " + filename);

        BufferedReader reader = new BufferedReader(new FileReader(filename));

        List<String> taxa = new ArrayList<>();

        String line = reader.readLine();
        while (line != null) {
            taxa.add(line);
            System.out.println("  " + line);
            line = reader.readLine();
        }
        return taxa;
    }

    private static String taxon(String line) {

        String[] pieces = line.split("taxon=\"");
        String taxon = pieces[1].substring(1, pieces[1].indexOf('"'));
        //taxon.replace('_', ' ');
        return taxon;
    }

    private static String newLowerUpper(String line, double lower, double upper) {
        String[] pieces = line.split("taxon=\"");
        String taxon = pieces[1].substring(1, pieces[1].indexOf('"'));

        StringBuilder builder = new StringBuilder();
        builder.append(pieces[0]);
        builder.append("taxon=\"@");
        builder.append(taxon);
        builder.append("\" lower=\"");
        builder.append(lower);
        builder.append("\" upper=\"");
        builder.append(upper);
        builder.append("\"/>");

        return builder.toString();
    }
}
