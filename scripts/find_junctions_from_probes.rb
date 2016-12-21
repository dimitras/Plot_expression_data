# USAGE:
# ruby find_junctions_from_probes.rb results/probes_mm9.bed.txt data/FINAL_master_list_of_junction_counts_MIN.lihong.txt results/selected_junctions_for_probes.txt

# Find the junctions that correspond to the probes.

require 'rubygems'
require 'csv'

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]

# read the list
probes = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1, {:col_sep => "\t"}) do |row|
	if row[0] != "" 
		start = row[1].to_i+1
		stop = row[2].to_i-1
		puts row[8]
		probes[row[8]] << [row[0], start, stop, row[6], row[8]]
	end
end

junctions = Hash.new { |h,k| h[k] = [] }
header = nil
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	if row[0] == "loc"
		header = row
	end
	if row[0] != "" && row[0] != "loc"
		chrom, start, stop = row[0].scan(/^junction:(chr\w{1,2}):(\d+)-(\d+)$/)
		if row[26].nil?
			gene = "NA"
		elsif row[26].include? ","
			gene = row[26].split(",")[0]
		else
			gene = row[26]
		end
		junctions[gene] << [row , chrom, start, stop, gene]
	end
end


selected_junctions = Hash.new { |h,k| h[k] = [] }
probes.each do |gene, list|
	puts gene
	puts list.inspect
	list.each do |row|
		junctions[gene].each do |erow|
			erow.flatten!(2)
			if row[0] == erow[27] # check chrom
				if (erow[28].to_i == row[1].to_i) && (erow[29].to_i == row[2].to_i)
					puts gene
					selected_junctions[gene] << erow.push(row[3]).push(row[1]).push(row[2])
				end
				# if (erow[28] <= row[2]) && (row[1] <= erow[29]) # check the overlap (exonstart<=probestop) && (probestart<=exonstop)
				# 	selected_junctions[gene] << erow.push(row[6]).push(row[1]).push(row[2])
				# end
			end
		end
	end
end


# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << header[0..25].flatten(2) + ["geneSymbol", "chrom", "start", "end", "probe", "start", "end"]
	selected_junctions.each do |gene, list|
		list.each do |row|
			csv << row.flatten(2)
		end
	end
end
