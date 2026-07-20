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

		fs_result = []
		indel_result = []
		ns_result = []

		if self.dna_xml_file:
			try:
				with open(self.dna_xml_file, 'r') as blastn_result_handle:
					blastn_records = NCBIXML.parse(blastn_result_handle)

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

									fs_dict_list = []
									pep_insert_dict_list = []
									pep_del_dict_list = []
									ns_dict_list = []

									if json_data[model_id]["model_param"].get("40494"):  # frameshifts
										for each_fs in list(json_data[model_id]["model_param"]["40494"]["param_value"].values()):
											original_aa, pos = self.parse_fsns(each_fs)
											fs_dict_list.append({
												"original_aa": original_aa,
												"aa_position": pos
											})

									if json_data[model_id]["model_param"].get("41344"):  # insertions into peptide seqs
										for eachpepin in list(json_data[model_id]["model_param"]["41344"]["param_value"].values()):
											if "_" in eachpepin:
												result = self.parse_indels(eachpepin)
												pep_insert_dict_list.append({
													"aa1": result[0],
													"pos1": int(result[1]),
													"aa2": result[2],
													"pos2": int(result[3]),
													"event": result[4],
													"deleted": result[5],
													"full_indel": eachpepin
													})
											else:
												result = self.parse_indels(eachpepin)
												pep_insert_dict_list.append({
													"aa1": result[0],
													"pos1": int(result[1]),
													"event": result[2],
													"aa2": result[3],
													"pos2": result[4],
													"deleted": result[5],
													"full_indel": eachpepin
													})

									if json_data[model_id]["model_param"].get("41342"): # deletions into peptide seqs
										for eachpepdel in list(json_data[model_id]["model_param"]["41342"]["param_value"].values()):
											if "_" in eachpepdel:
												result = self.parse_indels(eachpepdel)
												pep_del_dict_list.append({
													"aa1": result[0],
													"pos1": int(result[1]),
													"aa2": result[2],
													"pos2": int(result[3]),
													"event": result[4],
													"deleted": result[5],
													"full_indel": eachpepdel
													})
											else:
												result = self.parse_indels(eachpepdel)
												pep_del_dict_list.append({
													"aa1": result[0],
													"pos1": int(result[1]),
													"event": result[2],
													"aa2": result[3],
													"pos2": result[4],
													"deleted" : result[5],
													"full_indel": eachpepdel
													})
												
									if json_data[model_id]["model_param"].get("40394"):  # nonsense
										for eachns in list(json_data[model_id]["model_param"]["40394"]["param_value"].values()):
											original_aa, pos = self.parse_fsns(eachns)
											ns_dict_list.append({
												"original_aa": original_aa,
												"aa_position": pos
											})

									for hsp in alignment.hsps:
										query_seq =  hsp.query.replace('-', '')
										real_query_length = len(query_seq)
										sbjct_seq = hsp.sbjct.replace('-', '')
										real_sbjct_length = len(sbjct_seq)

										card_dna_ref = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]

										if fs_dict_list:
											fs_out = self.frameshift(hsp.query, hsp.sbjct, card_dna_ref, bnquery_def, param_type=json_data[model_id]["model_param"]["40494"]["param_type"], fs_dict_list=fs_dict_list)
										else:
											fs_out = None
										if pep_insert_dict_list or pep_del_dict_list:
											indel_out = self.indel(hsp.query, hsp.sbjct, card_dna_ref, bnquery_def, curated_in_list=pep_insert_dict_list, curated_del_list=pep_del_dict_list)
										else:
											indel_out = None
										if ns_dict_list:
											ns_out = self.nonsense(hsp.query, hsp.sbjct, card_dna_ref, bnquery_def, param_type=json_data[model_id]["model_param"]["40394"]["param_type"], ns_dict_list=ns_dict_list)
										else:
											ns_out = None

										# fetch mutations from MM
										if fs_out is not None:
											fs_result.append(fs_out)
										if indel_out is not None:
											indel_result.append(indel_out)
										if ns_out is not None:
											ns_result.append(ns_out)
								else:
									pass
						else:
							pass
			except FileNotFoundError as e:
				traceback.print_exc()
				logger.info("Skipping PVM extended mutation search...")
		else:
			logger.info("Skipping PVM extended mutation search...")

		with open(self.xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			for blast_record in blast_records:
				perfect = {}
				strict = {}
				loose = {}

				## filter MM results to only entries matching this blast_record's query
				bpquery_def = blast_record.query
				mutation_result = (fs_result or []) + (indel_result or []) + (ns_result or [])

				mutation_result_filtered = [
					m for m in mutation_result
					if m["query_def"].split()[0] in bpquery_def] if mutation_result else None
										
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

							try:
								for srv_result in self.single_resistance_variant(
									"PVM", snp_dict_list, hsp.query, hsp.sbjct_start, hsp.sbjct, orf_info, bpquery_def, 
									pred_genes_dict_prot=predicted_genes_dict_protein, sub_prot_dict=submitted_proteins_dict, real_sbjct_length=real_sbjct_length
									):
									mm_output = self.consolidate_mutations(
										self.input_type, 
										hit_id.decode(), 
										model_type="PVM", 
										srv=srv_result, 
										other_mutations=mutation_result_filtered, 
										hsp_bitscore=hsp.bits, 
										pass_val=true_pass_evalue)
									
									if not mm_output:
										continue
									
									mm_record = mm_output[0]
									curated_mutations = mm_record.get("curated_mutations", None)
									de_novo_mutations = mm_record.get("de_novo_mutations", None)

									if float(hsp.bits) >= float(true_pass_evalue):
										""" Strict hits """
										sinsidedict = {}
										sinsidedict["type_match"] = "Strict"
										if "eachs" in mm_record:
											sinsidedict["snp"] = mm_record["eachs"]
											sinsidedict["ast_source"] = self.get_ast_source(
												json_data[model_id], mm_record["eachs"])
											if "query_snps" in mm_record:
												sinsidedict["query_snp"] = mm_record["query_snps"]
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
											
										if curated_mutations is not None:
											sinsidedict["curated_mutations"] = '; '.join(', '.join(mutations) for mutations in curated_mutations.values())
											sinsidedict["curated_mutation_types"] = '; '.join(curated_mutations.keys())
										else:
											sinsidedict["curated_mutations"] = "n/a"
											sinsidedict["curated_mutation_types"] = "n/a"
										if de_novo_mutations is not None:
											sinsidedict["de_novo_mutations"] = '; '.join(', '.join(mutations) for mutations in de_novo_mutations.values())
											sinsidedict["de_novo_mutation_types"] = '; '.join(de_novo_mutations.keys())
										else:
											sinsidedict["de_novo_mutations"] = "n/a"
											sinsidedict["de_novo_mutation_types"] = "n/a"

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
										sinsidedict["sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
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
												if "orf_protein_sequence" in mm_record:
													sinsidedict["orf_prot_sequence"] = mm_record["orf_protein_sequence"]
												else:
													sinsidedict["orf_prot_sequence"] = "n/a"

											else:
												sinsidedict["orf_dna_sequence"] = ""
												sinsidedict["orf_prot_sequence"] = ""


										elif self.input_type == 'protein':
											sinsidedict["query_start"] = hsp.query_start
											sinsidedict["query_end"] = hsp.query_start + real_query_length
											sinsidedict["query_from"] = blast_record.query
											if "orf_protein_sequence" in mm_record:
												sinsidedict["orf_prot_sequence"] = mm_record["orf_protein_sequence"]
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
										if "eachs" in mm_record:
											slinsidedict["snp"] = mm_record["eachs"]
											slinsidedict["ast_source"] = self.get_ast_source(
												json_data[model_id], mm_record["eachs"])
											if "query_snps" in mm_record:
												slinsidedict["query_snp"] = mm_record["query_snps"]
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
											
										if curated_mutations is not None:
											slinsidedict["curated_mutations"] = '; '.join(', '.join(mutations) for mutations in curated_mutations.values())
											slinsidedict["curated_mutation_types"] = '; '.join(curated_mutations.keys())
										else:
											slinsidedict["curated_mutations"] = "n/a"
											slinsidedict["curated_mutation_types"] = "n/a"
										if de_novo_mutations is not None:
											slinsidedict["de_novo_mutations"] = '; '.join(', '.join(mutations) for mutations in de_novo_mutations.values())
											slinsidedict["de_novo_mutation_types"] = '; '.join(de_novo_mutations.keys())
										else:
											slinsidedict["de_novo_mutations"] = "n/a"
											slinsidedict["de_novo_mutation_types"] = "n/a"
											
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
												if "orf_protein_sequence" in mm_record:
													slinsidedict["orf_prot_sequence"] = mm_record["orf_protein_sequence"]
												else:
													slinsidedict["orf_prot_sequence"] = "n/a"

											else:
												slinsidedict["orf_dna_sequence"] = ""
												slinsidedict["orf_prot_sequence"] = ""

										elif self.input_type == 'protein':
											slinsidedict["query_start"] = hsp.query_start
											slinsidedict["query_end"] = hsp.query_start + real_query_length
											slinsidedict["query_from"] = blast_record.query
											if "orf_protein_sequence" in mm_record:
												slinsidedict["orf_prot_sequence"] = mm_record["orf_protein_sequence"]
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
								traceback.print_exc()
								logger.warning("Exception : {} -> {} -> Model({})".format(type(e), e, model_id))
								logger.warning("{} ---> hsp.bits: {} {} ? {}".format(json_data[model_id]["model_name"],hsp.bits,type(hsp.bits), type(true_pass_evalue)))

				blastResults = self.results(
					blastResults, blast_record.query, perfect, strict , loose, self.include_nudge)

			return blastResults
