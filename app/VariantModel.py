from app.Mutations import MutationsModule
from app.settings import *
from Bio.Blast import NCBIXML

class Variant(MutationsModule):
	"""Class for protein variant searches."""
	def __init__(self, input_type, loose, input_sequence, xml_file, dna_xml_file, working_directory, local_database=False, include_nudge=False):
		self.input_type = input_type
		self.loose = loose
		self.input_sequence = input_sequence
		self.xml_file = xml_file
		self.dna_xml_file = dna_xml_file
		self.output = {}
		self.working_directory = working_directory

		self.local_database = local_database
		self.data = data_path

		self.include_nudge = include_nudge

		if self.local_database:
			self.db = LOCAL_DATABASE
			self.data = LOCAL_DATABASE

	def __repr__(self):
		"""Returns Variant class full object."""
		return "Variant({}".format(self.__dict__)

	def run(self):
		blastResults = {}
		predicted_genes_dict = {}
		predicted_genes_dict_protein = {}
		submitted_proteins_dict = {}
		orf = 0

		if self.input_type == "contig":
			predicted_genes_dict = self.get_orf_dna_sequence(
				self.input_sequence, self.input_type)
			predicted_genes_dict_protein = self.get_orf_protein_sequence(
				self.input_sequence, self.input_type)
			
		# print("predicted proteins") # translated from the predicted GENES
		# print(predicted_genes_dict_protein)
		# print("=========================")
		# print("predicted genes")
		# print(predicted_genes_dict)
		# print("=========================")

		if self.input_type == "protein":
			submitted_proteins_dict = (
				self.get_submitted_protein_sequence(self.input_sequence))
		
		# print("submitted proteins (no predictions needed, just BLAST/DIAMOND)")
		# print(submitted_proteins_dict)
		# print("=========================")

		with open(os.path.join(self.data,"card.json")) as json_file:
			json_data = json.load(json_file)

		# if self.dna_xml_file == None:
		# 	pass
		# else:
		try:
			with open(self.dna_xml_file, 'r') as blastn_result_handle:
				blastn_records = NCBIXML.parse(blastn_result_handle)	
				for blastn_record in blastn_records:
					for alignment in blastn_record.alignments:
						align_title = alignment.title

						model_type_id = self.extract_nth_bar(align_title, 0)
						# logger.info("model_type_id: {} ".format(model_type_id))
						space_pos = align_title.index(' ')
						hit_id = align_title[0:space_pos]
						hit_id = hit_id.encode('ascii','replace')
						# print("debug 2:", hit_id)

						model_descrpt = align_title[align_title.index(' ')+1:]
						underscore_in_MD = model_descrpt.index('_')
						model_id = model_descrpt[0:underscore_in_MD]
						seq_in_model = model_descrpt[underscore_in_MD+1: model_descrpt.index(' ')]
						pass_value = self.extract_nth_bar(alignment.title, 1)
						# logger.info("pass_value: {}".format(pass_value))

						if model_type_id == 40293:
							try:
								true_pass_evalue = float(pass_value)
							except ValueError:
								true_pass_evalue = float(
									pass_value[0:pass_value.find(' ')])

							logger.info("mutation | model_type_id = " + str(align_title))
							init = 0
							fsl = []
							temp = ""
							
							evalue_fs = self.extract_nth_bar(align_title, 2)
							fsl = evalue_fs.split(',')
							
							for hsp in alignment.hsps:
								query_seq =  hsp.query.replace('-', '')
								real_query_length = len(query_seq)
								sbjct_seq = hsp.sbjct.replace('-', '')
								real_sbjct_length = len(sbjct_seq)

								card_dna_ref = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]

								if "Frameshift" in align_title:	
									self.frameshift(fsl, hsp.query, hsp.sbjct, card_dna_ref)
		except FileNotFoundError as e:
			logger.info("Skipping frameshift search...")
			pass

		with open(self.xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			for blast_record in blast_records:
				perfect = {}
				strict = {}
				loose = {}
				for alignment in blast_record.alignments:
					align_title = alignment.title
					orf_info = blast_record.query.encode('ascii','replace')
					c = 0
					barc = 0
					for eachc in orf_info:
						if barc >= 6:
							break
						elif eachc == '|':
							barc += 1
							c += 1
						else:
							c += 1
					orf_from = orf_info[c:]

					model_type_id = self.extract_nth_bar(align_title, 0)
					# print(model_type_id)
					# logger.info("model_type_id: {} ".format(model_type_id))
					space_pos = align_title.index(' ')
					hit_id = align_title[0:space_pos]
					hit_id = hit_id.encode('ascii','replace')
					# print("debug 2:", hit_id)

					model_descrpt = align_title[align_title.index(' ')+1:]
					underscore_in_MD = model_descrpt.index('_')
					model_id = model_descrpt[0:underscore_in_MD]
					seq_in_model = model_descrpt[underscore_in_MD+1: model_descrpt.index(' ')]
					pass_value = self.extract_nth_bar(alignment.title, 1)
					# logger.info("pass_value: {}".format(pass_value))

					if model_type_id == 40293:
						try:
							true_pass_evalue = float(pass_value)
						except ValueError:
							true_pass_evalue = float(
								pass_value[0:pass_value.find(' ')])

						logger.info("mutation | model_type_id = " + str(align_title))
						init = 0
						snpl = []
						temp = ""

						evalue_snp = self.extract_nth_bar(align_title, 2)
						# evalue_snp_dec = evalue_snp
						snpl = evalue_snp.split(',')
						
						# print("\nblastp align title:", align_title)
						# print(snpl)
						# print()

						for hsp in alignment.hsps:
							query_seq =  hsp.query.replace('-', '')
							real_query_length = len(query_seq)
							sbjct_seq = hsp.sbjct.replace('-', '')
							real_sbjct_length = len(sbjct_seq)

							# print(orf_info,"\n",hsp)

							# print(self.extract_nth_bar(orf_info.decode(), 0))
							# print(self.extract_nth_bar(orf_info.decode(), 1))
							# print(self.extract_nth_bar(orf_info.decode(), 2))
							# print(orf_from.decode())

							self.single_resistance_variant(
								predicted_genes_dict_protein, submitted_proteins_dict, snpl, real_sbjct_length, hsp.query, hsp.sbjct_start, hsp.sbjct, orf_info
								)
							
							mm_output = self.consolidate_mutations()
							# print(mm_output)
							# for s in mm_output:
							# 	print(s)
							
							if mm_output != None:
								for loaded_snp in mm_output:
									if len(loaded_snp) > 0: ## to weed out SNPs without ORFs or qry entries
										# print("loaded snp:", loaded_snp)
										try:
											# print("debug 2:", hit_id, ",", hsp.query[loaded_snp["qry"]], ",", loaded_snp["qry"], ",", loaded_snp["chan"])
											if float(hsp.bits) >= float(true_pass_evalue): # if the hit passes its bitscore cut off (but isn't perfect)
												# print("hsp bits:", hsp.bits)
												sinsidedict = {}
												sinsidedict["type_match"] = "Strict"
												if "eachs" in loaded_snp:
													sinsidedict["snp"] = loaded_snp["eachs"]
													sinsidedict["query_snp"] = loaded_snp["query_snps"]
													sinsidedict["orf_strand"] = self.extract_nth_bar(orf_info.decode(), 0)
													sinsidedict["orf_start"] = self.extract_nth_bar(orf_info.decode(), 1)
													sinsidedict["orf_end"] = self.extract_nth_bar(orf_info.decode(), 2)
													sinsidedict["orf_from"] = orf_from.decode()
												else:
													sinsidedict["snp"] = "n/a"
													sinsidedict["query_snp"] = "n/a"
													sinsidedict["orf_strand"] = "n/a"
													sinsidedict["orf_start"] = "n/a"
													sinsidedict["orf_end"] = "n/a"
													sinsidedict["orf_from"] = "n/a"
												if "curated_fs" in loaded_snp:
													sinsidedict["curated_fs"] = loaded_snp["curated_fs"]
												else:
													sinsidedict["curated_fs"] = "n/a"
												if "denovo_fs" in loaded_snp:
													sinsidedict["denovo_fs"] = loaded_snp["denovo_fs"]
												else:
													sinsidedict["denovo_fs"] = "n/a"
												sinsidedict["model_name"] = json_data[model_id]["model_name"]
												sinsidedict["model_type"] = json_data[model_id]["model_type"]
												sinsidedict["model_type_id"] = model_type_id
												sinsidedict["model_id"] = model_id
												sinsidedict["pass_evalue"] = "n/a"
												sinsidedict["pass_bitscore"] = pass_value
												sinsidedict["ARO_accession"] = json_data[model_id]["ARO_accession"]
												sinsidedict["ARO_name"] = json_data[model_id]["ARO_name"]
												sinsidedict["ARO_category"] = json_data[model_id]["ARO_category"]
												sinsidedict["evalue"] = hsp.expect
												sinsidedict["max_identities"] = hsp.identities
												sinsidedict["bit_score"] = hsp.bits
												sinsidedict["cvterm_id"]  = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
												sinsidedict["query"] = hsp.query
												sinsidedict["match"] = hsp.match
												sinsidedict["sequence_from_db"] = hsp.sbjct
												sinsidedict["sequence_from_broadstreet"]	= json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
												sinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
												if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
													sinsidedict["partial"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["partial"]
												else:
													sinsidedict["partial"] = "0"

												if self.input_type == 'contig':
													sinsidedict["query_start"] = self.extract_nth_hash(
														orf_info.decode(), 1) + (hsp.query_start - 1)*3
													sinsidedict["query_end"] = self.extract_nth_hash(
														orf_info.decode(), 1) + (hsp.query_start - 1)*3 + real_query_length*3 - 1
													sinsidedict["orf_strand"] = self.extract_nth_hash(
														orf_info.decode(), 3)
													sinsidedict["orf_start"] = self.extract_nth_hash(
														orf_info.decode(), 1)
													sinsidedict["orf_end"] = self.extract_nth_hash(
														orf_info.decode(), 2)
													sinsidedict["orf_from"] = self.extract_nth_hash(
														orf_info.decode(), 0)
													sinsidedict["hit_start"] = (
														hsp.sbjct_start-1)*3
													sinsidedict["hit_end"] = (
														hsp.sbjct_end)*3

													if orf_info.decode().split(' # ')[0] in predicted_genes_dict:
														sinsidedict["orf_dna_sequence"] = predicted_genes_dict[orf_info.decode().split(' # ')[0]]
														# sinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orf_info.decode().split(' # ')[0]], generic_dna).translate(table=11)).strip("*")
														if "orf_protein_sequence" in loaded_snp:
															sinsidedict["orf_prot_sequence"] = loaded_snp["orf_protein_sequence"]
														else:
															sinsidedict["orf_prot_sequence"] = "n/a"
														# print("debug 1:", value)

													else:
														sinsidedict["orf_dna_sequence"] = ""
														sinsidedict["orf_prot_sequence"] = ""


												elif self.input_type == 'protein':
													sinsidedict["query_start"] = hsp.query_start
													sinsidedict["query_end"] = hsp.query_start + real_query_length
													sinsidedict["query_from"] = blast_record.query
													if "orf_protein_sequence" in loaded_snp:
														sinsidedict["orf_prot_sequence"] = loaded_snp["orf_protein_sequence"]
													else:
														sinsidedict["orf_prot_sequence"] = "n/a"
													# print("debug 2:", value)

													sinsidedict["hit_start"] = ""
													sinsidedict["hit_end"] = ""

												elif self.input_type == 'read':
													pass

												sinsidedict["perc_identity"] = float(format(
													float(sinsidedict["max_identities"]*100) / len(sinsidedict["query"]), '.2f'))

												strict["{}|hsp_num:{}".format(
													hit_id.decode(), init)] = sinsidedict
												init += 1

											else:
												# print("hsp bits:", hsp.bits)
												slinsidedict = {}
												slinsidedict["type_match"] = "Loose"
												if "eachs" in loaded_snp:
													slinsidedict["snp"] = loaded_snp["eachs"]
													slinsidedict["query_snp"] = loaded_snp["query_snps"]
													slinsidedict["orf_strand"] = self.extract_nth_bar(orf_info.decode(), 0)
													slinsidedict["orf_start"] = self.extract_nth_bar(orf_info.decode(), 1)
													slinsidedict["orf_end"] = self.extract_nth_bar(orf_info.decode(), 2)
													slinsidedict["orf_from"] = orf_from.decode()
												else:
													slinsidedict["snp"] = "n/a"
													slinsidedict["query_snp"] = "n/a"
													slinsidedict["orf_strand"] = "n/a"
													slinsidedict["orf_start"] = "n/a"
													slinsidedict["orf_end"] = "n/a"
													slinsidedict["orf_from"] = "n/a"
												if "curated_fs" in loaded_snp:
													slinsidedict["curated_fs"] = loaded_snp["curated_fs"]
												else:
													slinsidedict["curated_fs"] = "n/a"
												if "denovo_fs" in loaded_snp:
													slinsidedict["denovo_fs"] = loaded_snp["denovo_fs"]
												else:
													slinsidedict["denovo_fs"] = "n/a"
												slinsidedict["model_name"] = json_data[model_id]["model_name"]
												slinsidedict["model_type"] = json_data[model_id]["model_type"]
												slinsidedict["model_type_id"] = model_type_id
												slinsidedict["pass_evalue"] = "n/a"
												slinsidedict["pass_bitscore"] = pass_value
												slinsidedict["model_id"] = model_id
												slinsidedict["ARO_accession"] = json_data[model_id]["ARO_accession"]
												slinsidedict["ARO_name"] = json_data[model_id]["ARO_name"]
												slinsidedict["ARO_category"] = json_data[model_id]["ARO_category"]
												slinsidedict["evalue"] = hsp.expect
												slinsidedict["bit_score"] = hsp.bits
												slinsidedict["max_identities"] = hsp.identities
												slinsidedict["cvterm_id"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
												slinsidedict["query"] = hsp.query
												slinsidedict["match"] = hsp.match
												slinsidedict["sequence_from_db"] = hsp.sbjct
												slinsidedict["sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
												slinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
												if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
													slinsidedict["partial"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["partial"]
												else:
													slinsidedict["partial"] = "0"

												if self.input_type == 'contig':
													slinsidedict["query_start"] = self.extract_nth_hash(
														orf_info.decode(), 1) + (hsp.query_start - 1)*3
													slinsidedict["query_end"] = self.extract_nth_hash(
														orf_info.decode(), 1) + (hsp.query_start - 1)*3 + real_query_length*3 - 1
													slinsidedict["orf_strand"] = self.extract_nth_hash(
														orf_info.decode(), 3)
													slinsidedict["orf_start"] = self.extract_nth_hash(
														orf_info.decode(), 1)
													slinsidedict["orf_end"] = self.extract_nth_hash(
														orf_info.decode(), 2)
													slinsidedict["orf_from"] = self.extract_nth_hash(
														orf_info.decode(), 0)
													slinsidedict["hit_start"] = (
														hsp.sbjct_start-1)*3
													slinsidedict["hit_end"] = (
														hsp.sbjct_end)*3

													if orf_info.decode().split(' # ')[0] in predicted_genes_dict:
														slinsidedict["orf_dna_sequence"] = predicted_genes_dict[orf_info.decode().split(' # ')[0]]
														# slinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orf_info.decode().split(' # ')[0]], generic_dna).translate(table=11)).strip("*")
														if "orf_protein_sequence" in loaded_snp:
															slinsidedict["orf_prot_sequence"] = loaded_snp["orf_protein_sequence"]
														else:
															slinsidedict["orf_prot_sequence"] = "n/a"
														# print("debug 3:", value)

													else:
														slinsidedict["orf_dna_sequence"] = ""
														slinsidedict["orf_prot_sequence"] = ""

												elif self.input_type == 'protein':
													slinsidedict["query_start"] = hsp.query_start
													slinsidedict["query_end"] = hsp.query_start + real_query_length
													slinsidedict["query_from"] = blast_record.query
													if "orf_protein_sequence" in loaded_snp:
														slinsidedict["orf_prot_sequence"] = loaded_snp["orf_protein_sequence"]
													else:
														slinsidedict["orf_prot_sequence"] = "n/a"
													# print("debug 4:", value)

													# print("debug:", slinsidedict["query_start"], slinsidedict["query_end"], slinsidedict["query_from"])
													slinsidedict["hit_start"] = ""
													slinsidedict["hit_end"] = ""

												elif self.input_type == 'read':
													pass

												slinsidedict["perc_identity"] = float(format(
													float(slinsidedict["max_identities"]*100) / len(slinsidedict["query"]), '.2f'))
												loose["{}|hsp_num:{}".format(
													hit_id.decode(), init)] = slinsidedict

												init += 1
											# print(slinsidedict)

										except Exception as e:
											logger.warning("Exception : {} -> {} -> Model({})".format(type(e), e, model_id))
											logger.warning("{} ---> hsp.bits: {} {} ? {}".format(json_data[model_id]["model_name"],hsp.bits,type(hsp.bits), type(true_pass_evalue)))
									# else:
									# 	print(hsp.query[loaded_snp["qry"]], "+", loaded_snp["chan"])

				blastResults = self.results(blastResults, blast_record.query, perfect, strict , loose, self.include_nudge)
				# print(blastResults)

		return blastResults