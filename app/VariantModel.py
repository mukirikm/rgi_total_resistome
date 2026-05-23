from app.Mutations import MutationsModule
from app.settings import *
from Bio.Blast import NCBIXML
import traceback

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

		if self.input_type == "protein":
			submitted_proteins_dict = (
				self.get_submitted_protein_sequence(self.input_sequence))

		with open(os.path.join(self.data,"card.json")) as json_file:
			json_data = json.load(json_file)

		try:
			with open(self.dna_xml_file, 'r') as blastn_result_handle:
				blastn_records = NCBIXML.parse(blastn_result_handle)
				fs_result = []
				indel_result = []

				for blastn_record in blastn_records:
					bnquery_def = blastn_record.query
					if blastn_record.alignments:
						for alignment in blastn_record.alignments:	
							align_title = alignment.title
							model_type_id = self.extract_nth_bar(align_title, 0)
							space_pos = align_title.index(' ')
							hit_id = align_title[0:space_pos]
							hit_id = hit_id.encode('ascii','replace')
							model_descrpt = align_title[align_title.index(' ')+1:]
							underscore_in_MD = model_descrpt.index('_')
							model_id = model_descrpt[0:underscore_in_MD]
							seq_in_model = model_descrpt[underscore_in_MD+1: model_descrpt.index(' ')]
							pass_value = self.extract_nth_bar(alignment.title, 1)
							
							if model_type_id == 40293:
								try:
									true_pass_evalue = float(pass_value)
								except ValueError:
									true_pass_evalue = float(
										pass_value[0:pass_value.find(' ')])

								fsl = []
								fs_dict_list = []
								
								evalue_fs = self.extract_nth_bar(align_title, 2)
								fsl = evalue_fs.split(',')

								## grabbing curated frameshifts from blast XML (change to CARD JSON as input later?)
								for each_fs in fsl:
									position = int(
										''.join(filter(str.isdigit, each_fs)))

									original = (each_fs.split(
										''.join(filter(str.isdigit, each_fs))))
									
									fs_dict_list.append(
										{"original_aa": original[0], "aa_position": position})
								
								for hsp in alignment.hsps:
									query_seq =  hsp.query.replace('-', '')
									real_query_length = len(query_seq)
									sbjct_seq = hsp.sbjct.replace('-', '')
									real_sbjct_length = len(sbjct_seq)

									card_dna_ref = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]

									# fetch mutations from MM
									fs_result.append(self.frameshift(hsp.query, hsp.sbjct, card_dna_ref, bnquery_def, fs_dict_list=fs_dict_list))
									
									indel_result.append(self.indel(hsp.query, hsp.sbjct, card_dna_ref, bnquery_def))
					else:
						fs_result.append({"query_def": bnquery_def, "has_fs": False})
						indel_result.append({"query_def": bnquery_def, "has_indel": False})
				
		except FileNotFoundError as e:
			fs_result = None
			logger.info("Skipping PVM frameshift search...")
			pass

		with open(self.xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			for blast_record in blast_records:
				perfect = {}
				strict = {}
				loose = {}

				## filter fs_result to only entries matching this blast_record's query
				bpquery_def = blast_record.query
				fs_result_filtered = [
					f for f in fs_result 
					if f["query_def"].split()[0] in bpquery_def] if fs_result else None
		
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
					space_pos = align_title.index(' ')
					hit_id = align_title[0:space_pos]
					hit_id = hit_id.encode('ascii', 'replace')
					model_descrpt = align_title[align_title.index(' ')+1:]
					underscore_in_MD = model_descrpt.index('_')
					model_id = model_descrpt[0:underscore_in_MD]
					seq_in_model = model_descrpt[underscore_in_MD +
                                                 1: model_descrpt.index(' ')]
					pass_value = self.extract_nth_bar(alignment.title, 1)
					model_type_id = self.extract_nth_bar(align_title, 0)
					# logger.info("model_type_id: {} ".format(model_type_id))
					space_pos = align_title.index(' ')
					hit_id = align_title[0:space_pos]
					hit_id = hit_id.encode('ascii','replace')

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

						# logger.info("mutation | model_type_id = " + str(align_title))
						init = 0
						snpl = []
						snp_dict_list = []
						temp = ""

						evalue_snp = self.extract_nth_bar(align_title, 2)
						# evalue_snp_dec = evalue_snp
						snpl = evalue_snp.split(',')

						for each_snp in snpl:
							position = int(
                                ''.join(filter(str.isdigit, each_snp)))
							
							original_change = (each_snp.split(
                                ''.join(filter(str.isdigit, each_snp))))
							
							snp_dict_list.append(
                                {"original": original_change[0], "change": original_change[-1], "position": position})

						for hsp in alignment.hsps:
							query_seq =  hsp.query.replace('-', '')
							real_query_length = len(query_seq)
							sbjct_seq = hsp.sbjct.replace('-', '')
							real_sbjct_length = len(sbjct_seq)

							for srv_result in self.single_resistance_variant(
								"PVM", snp_dict_list, hsp.query, hsp.sbjct_start, hsp.sbjct, orf_info, bpquery_def, 
								pred_genes_dict_prot=predicted_genes_dict_protein, sub_prot_dict=submitted_proteins_dict, real_sbjct_length=real_sbjct_length
								):
								mm_output = self.consolidate_mutations(
									self.input_type, 
									hit_id.decode(), 
									model_type="pvm", 
									srv=srv_result, 
									fs=fs_result_filtered, 
									hsp_bitscore=hsp.bits, 
									pass_val=true_pass_evalue)

								try:
									if mm_output:
										for loaded_snp in mm_output:
											if float(hsp.bits) >= float(true_pass_evalue):
												""" Strict hits """
												sinsidedict = {}
												sinsidedict["type_match"] = "Strict"
												if "eachs" in loaded_snp:
													sinsidedict["snp"] = loaded_snp["eachs"]
													sinsidedict["ast_source"] = self.get_ast_source(
														json_data[model_id], loaded_snp["eachs"])
													if "query_snps" in loaded_snp:
														sinsidedict["query_snp"] = loaded_snp["query_snps"]
													else:
														sinsidedict["query_snp"] = "n/a"
													sinsidedict["orf_strand"] = self.extract_nth_bar(orf_info.decode(), 0)
													sinsidedict["orf_start"] = self.extract_nth_bar(orf_info.decode(), 1)
													sinsidedict["orf_end"] = self.extract_nth_bar(orf_info.decode(), 2)
													sinsidedict["orf_from"] = self.trim_after_last_underscore(orf_from.decode())
												else:
													sinsidedict["snp"] = "n/a"
													sinsidedict["ast_source"] = "n/a"
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
													sinsidedict["orf_from"] = self.trim_after_last_underscore(self.extract_nth_hash(
														orf_info.decode(), 0))
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
												""" Loose hits """
												slinsidedict = {}
												slinsidedict["type_match"] = "Loose"
												if "eachs" in loaded_snp:
													slinsidedict["snp"] = loaded_snp["eachs"]
													slinsidedict["ast_source"] = self.get_ast_source(
														json_data[model_id], loaded_snp["eachs"])
													if "query_snps" in loaded_snp:
														slinsidedict["query_snp"] = loaded_snp["query_snps"]
													else:
														slinsidedict["query_snp"] = "n/a"
													slinsidedict["orf_strand"] = self.extract_nth_bar(orf_info.decode(), 0)
													slinsidedict["orf_start"] = self.extract_nth_bar(orf_info.decode(), 1)
													slinsidedict["orf_end"] = self.extract_nth_bar(orf_info.decode(), 2)
													slinsidedict["orf_from"] = self.trim_after_last_underscore(orf_from.decode())
												else:
													slinsidedict["snp"] = "n/a"
													slinsidedict["ast_source"] = "n/a"
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
													slinsidedict["orf_from"] = self.trim_after_last_underscore(self.extract_nth_hash(
														orf_info.decode(), 0))
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
													slinsidedict["hit_start"] = ""
													slinsidedict["hit_end"] = ""

												elif self.input_type == 'read':
													pass

												slinsidedict["perc_identity"] = float(format(
													float(slinsidedict["max_identities"]*100) / len(slinsidedict["query"]), '.2f'))
												loose["{}|hsp_num:{}".format(
													hit_id.decode(), init)] = slinsidedict

												init += 1

								except Exception as e:
									traceback.print_exc() # for karyn
									logger.warning("Exception : {} -> {} -> Model({})".format(type(e), e, model_id))
									logger.warning("{} ---> hsp.bits: {} {} ? {}".format(json_data[model_id]["model_name"],hsp.bits,type(hsp.bits), type(true_pass_evalue)))
								else:
									pass
				blastResults = self.results(
					blastResults, blast_record.query, perfect, strict , loose, self.include_nudge)

			return blastResults