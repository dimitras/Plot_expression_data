# USAGE:
# ruby find_uniongenes_in_allgenes.rb data/union3DEcomparisons_fc1.5pv0.1_union.csv data/FINAL_master_list_of_gene_counts_MIN.lihong_renamed.txt results/uniongenes_fullexpression.csv

# Filter the union genes in all genes list. 

require 'rubygems'
require 'csv'

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]

# read the list
uniongenes = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1) do |row|
	if row[0] != "" && row[0] != "id"
		uniongenes[row[0]] = row
	end
end

header = nil
allgenes = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	if row[0] == "id"
		header = row
	end
	if row[0] != "" && row[0] != "id"
		allgenes[row[0]] = row
	end
end


ugenes_wfullexpr = Hash.new { |h,k| h[k] = [] }
allgenes.each do |gene, row|
	if uniongenes.has_key?(gene)
		ugenes_wfullexpr[gene] = row
	end
end


# output
CSV.open(ofile, "wb") do |csv|
	csv << header[0..24]
	ugenes_wfullexpr.each do |gene, row|
		csv << row[0..24]
	end
end
