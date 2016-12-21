# USAGE:
# ruby find_exons_from_probes.rb probes_mm9.bed.txt results/filtered_exons_matching_probes.txt results/selected_exons_for_probes2.txt

# Find the exons that correspond to the probes.

require 'rubygems'
require 'csv'

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]

# read the list
exons = Hash.new { |h,k| h[k] = [] }
header = nil
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	if row[0] == "loci"
		header = row
	end
	if row[0] != "" && row[0] != "loci"
		gene = row[26].split(",")[0]
		exons[gene] << row
	end
end

probes = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1, {:col_sep => "\t"}) do |row|
	if row[0] != "" 
		probes[row[9]] << row
	end
end


selected_exons = Hash.new { |h,k| h[k] = [] }
probes.each do |gene, list|
	list.each do |row|
		exons[gene].each do |erow|
			if row[0] == erow[27] # check chrom
				if (erow[28] <= row[2]) && (row[1] <= erow[29]) # check the overlap (exonstart<=probestop) && (probestart<=exonstop)
					selected_exons[gene] << erow.push(row[8]).push(row[6])
				end
			end
		end
	end
end


# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << header[0..25].flatten(2) + ["geneSymbol", "chrom", "start", "end", "exon", "probe"]
	selected_exons.each do |gene, list|
		list.each do |row|
			csv << row.flatten(2)
		end
	end
end
