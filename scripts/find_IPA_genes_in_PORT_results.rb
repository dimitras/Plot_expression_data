# USAGE:
# ruby scripts/find_IPA_genes_in_PORT_results.rb IPA/pathways/ data/PORT/FINAL_master_list_of_gene_counts_MIN.lihong_renamed.txt results/Heatmaps/PORT.expressions.for.IPA/

# Find the pathway genes found from IPA, in the PORT results.

require 'rubygems'
require 'csv'

idir = ARGV[0]
ifile = ARGV[1]
ofile = ARGV[2]


Dir.foreach(idir) do |ifile1|
	next if !ifile1.include? ".csv"

	ofname = File.basename(ifile1, '.csv')

	genes = Hash.new { |h,k| h[k] = [] }
	CSV.foreach(idir+ifile1) do |row|
		if row[0] != "" && row[0] != "Symbol"
			genes[row[2]] = row[0]
		end
	end

	retrieved_genes = Hash.new { |h,k| h[k] = [] }
	header = nil
	CSV.foreach(ifile, {:col_sep => "\t"}) do |row|
		if row[0] == "id"
			header = row[0..25]
		end
		if row[0] != "" && row[0] != "id" && genes.has_key?(row[0])
			retrieved_genes[row[0]] = row[1..25]
		end
	end

		# output
	CSV.open(ofile+ofname+".fromPORT.txt", "wb", {:col_sep => "\t"}) do |csv|
		csv << [header, "EnsembleID"].flatten!(1)
		retrieved_genes.each do |gene, row|
			csv << [genes[gene], row, gene].flatten!(1)
		end
	end

end




