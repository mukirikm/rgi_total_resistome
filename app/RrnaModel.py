from app.Mutations import MutationsModule
from app.settings import *
from Bio.Blast import NCBIXML


class Rrna(MutationsModule):
    """Class for ribosomal RNA searches."""

    def __init__(self, input_file, output_file, db, xml, loose, local_database=False, include_nudge=False, num_threads=32):
        self.input_file = input_file
        self.output_file = output_file
        self.db = db
        self.xml_file = xml
        self.loose = loose

        self.local_database = local_database
        self.data = data_path

        self.include_nudge = include_nudge
        self.num_threads = num_threads

        if self.local_database:
            # self.db = LOCAL_DATABASE
            self.data = LOCAL_DATABASE

    def __repr__(self):
        """Returns Ribosomal RNA class full object."""
        return "Rrna({}".format(self.__dict__)

    def sequence_orientation(self, end, start):
        if end > start:
            return "+"
        else:
            return "-"

    def run(self):
        blastResults = {}
        with open(os.path.join(self.data, "card.json")) as json_file:
            json_data = json.load(json_file)

        with open(os.path.join(self.xml_file), 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)

            for blast_record in blast_records:
                perfect = {}
                strict = {}
                loose = {}

                ## filter fs_result to only entries matching this blast_record's query
                bpquery_def = blast_record.query
				# fs_result_filtered = [f for f in fs_result if f["query_def"].split()[0] in bpquery_def] if fs_result else None

                for alignment in blast_record.alignments:
                    align_title = alignment.title
                    orf_info = blast_record.query

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
                    orf_info_str = str(orf_info).split(" | ")
                    model_type_id = int(orf_info_str[1].split(":")[1].strip())
                    # logger.debug("model_type_id: {} ".format(model_type_id))

                    space_pos = align_title.index(' ')

                    hit_id = align_title[0:space_pos]
                    hit_id = hit_id.encode('ascii', 'replace')

                    model_descrpt = "?model_descrpt?"
                    model_info = orf_info_str[0].split("_")
                    model_id = model_info[0].strip()
                    seq_in_model = model_info[1].strip()
                    pass_value = orf_info_str[2].split(":")[1].strip()

                    # logger.debug("model_id: {}".format(model_id))
                    # logger.debug("pass_value: {}".format(pass_value))

                    if model_type_id == 40295:
                        predicted_genes_dict_protein = None
                        submitted_proteins_dict = None

                        true_pass_evalue = float(pass_value)

                        init = 0
                        evalue_snp = orf_info_str[3].split(":")[1].strip()
                        snpl = []
                        snp_dict_list = []
                        temp = ""
                        snpl = evalue_snp.split(',')

                        for each_snp in snpl:
                            snp_dict_list.append(
                                {"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})

                        for hsp in alignment.hsps:
                            query_seq = hsp.query.replace('-', '')
                            real_query_length = len(query_seq)
                            sbjct_seq = hsp.sbjct.replace('-', '')
                            real_sbjct_length = len(sbjct_seq)
                            strand = self.sequence_orientation(
                                hsp.sbjct_end, hsp.sbjct_start)
                            
                            for srv_result in self.single_resistance_variant(
								"RGV", snp_dict_list, hsp.query, hsp.sbjct_start, hsp.sbjct, orf_info, bpquery_def, 
                                hsp_query_start=hsp.query_start, hsp_query_end=hsp.query_end, real_qry_length=real_query_length, real_sbjct_length=real_sbjct_length, strand=strand
								):
                                    try:
                                        if float(hsp.bits) >= float(true_pass_evalue):
                                            sinsidedict = {}
                                            sinsidedict["type_match"] = "Strict"
                                            sinsidedict["ast_source"] = self.get_ast_source(
                                                json_data[model_id], srv_result["eachs"])
                                            sinsidedict["snp"] = srv_result["eachs"]
                                            sinsidedict["query_snp"] = srv_result["query_snps"]
                                            sinsidedict["orf_strand"] = strand
                                            sinsidedict["orf_start"] = hsp.sbjct_start
                                            sinsidedict["orf_end"] = hsp.sbjct_end

                                            sinsidedict["_orf_strand"] = self.extract_nth_bar(
                                                orf_info, 0)
                                            sinsidedict["_orf_start"] = self.extract_nth_bar(
                                                orf_info, 1)
                                            sinsidedict["_orf_end"] = self.extract_nth_bar(
                                                orf_info, 2)

                                            sinsidedict["orf_from"] = alignment.hit_def
                                            sinsidedict["strand"] = strand
                                            sinsidedict["hit_def"] = alignment.hit_def
                                            sinsidedict["sbjct_start"] = hsp.sbjct_start
                                            sinsidedict["sbjct_end"] = hsp.sbjct_end
                                            sinsidedict["query_start"] = hsp.query_start
                                            sinsidedict["query_end"] = hsp.query_end
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
                                            sinsidedict["cvterm_id"] = json_data[model_id]["model_sequences"][
                                                "sequence"][seq_in_model]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
                                            sinsidedict["query"] = hsp.query
                                            sinsidedict["match"] = hsp.match
                                            sinsidedict["sbjct"] = hsp.sbjct
                                            sinsidedict["sequence_from_db"] = ""
                                            sinsidedict["orf_dna_sequence"] = hsp.sbjct
                                            sinsidedict["orf_prot_sequence"] = ""
                                            sinsidedict["hit_start"] = (
                                                hsp.sbjct_start-1)*3
                                            sinsidedict["hit_end"] = (
                                                hsp.sbjct_end)*3

                                            sinsidedict["sequence_from_broadstreet"] = json_data[model_id][
                                                "model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
                                            sinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id][
                                                "model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
                                            if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
                                                sinsidedict["partial"] = json_data[model_id]["model_sequences"][
                                                    "sequence"][seq_in_model]["dna_sequence"]["partial"]
                                            else:
                                                sinsidedict["partial"] = "0"
                                            sinsidedict["perc_identity"] = float(format(
                                                float(sinsidedict["max_identities"]*100) / len(sinsidedict["query"]), '.2f'))

                                            strict["{}|hsp_num:{}".format(
                                                hit_id.decode(), init)] = sinsidedict
                                            init += 1

                                        else:
                                            slinsidedict = {}
                                            slinsidedict["type_match"] = "Loose"
                                            slinsidedict["ast_source"] = self.get_ast_source(
                                                json_data[model_id], srv_result["eachs"])
                                            slinsidedict["snp"] = srv_result["eachs"]
                                            slinsidedict["query_snp"] = srv_result["query_snps"]
                                            slinsidedict["orf_strand"] = strand
                                            slinsidedict["orf_start"] = hsp.sbjct_start
                                            slinsidedict["orf_end"] = hsp.sbjct_end

                                            slinsidedict["_orf_strand"] = self.extract_nth_bar(
                                                orf_info, 0)
                                            slinsidedict["_orf_start"] = self.extract_nth_bar(
                                                orf_info, 1)
                                            slinsidedict["_orf_end"] = self.extract_nth_bar(
                                                orf_info, 2)

                                            slinsidedict["orf_from"] = alignment.hit_def
                                            slinsidedict["strand"] = strand
                                            slinsidedict["hit_def"] = alignment.hit_def
                                            slinsidedict["sbjct_start"] = hsp.sbjct_start
                                            slinsidedict["sbjct_end"] = hsp.sbjct_end
                                            slinsidedict["query_start"] = hsp.query_start
                                            slinsidedict["query_end"] = hsp.query_end
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
                                            slinsidedict["cvterm_id"] = json_data[model_id]["model_sequences"][
                                                "sequence"][seq_in_model]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
                                            slinsidedict["query"] = hsp.query
                                            slinsidedict["match"] = hsp.match
                                            slinsidedict["sequence_from_db"] = ""
                                            slinsidedict["sbjct"] = hsp.sbjct
                                            slinsidedict["orf_dna_sequence"] = hsp.sbjct
                                            slinsidedict["orf_prot_sequence"] = ""
                                            slinsidedict["hit_start"] = (
                                                hsp.sbjct_start-1)*3
                                            slinsidedict["hit_end"] = (
                                                hsp.sbjct_end)*3

                                            slinsidedict["sequence_from_broadstreet"] = json_data[model_id][
                                                "model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
                                            slinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id][
                                                "model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
                                            if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
                                                slinsidedict["partial"] = json_data[model_id]["model_sequences"][
                                                    "sequence"][seq_in_model]["dna_sequence"]["partial"]
                                            else:
                                                slinsidedict["partial"] = "0"
                                            slinsidedict["perc_identity"] = float(format(
                                                float(slinsidedict["max_identities"]*100) / len(slinsidedict["query"]), '.2f'))

                                            loose["{}|hsp_num:{}".format(
                                                hit_id.decode(), init)] = slinsidedict
                                            init += 1
                                    except Exception as e:
                                        logger.warning(
                                            "Exception : {} -> {} -> Model({})".format(type(e), e, model_id))
                                        logger.warning("{} ---> hsp.bits: {} {} ? {}".format(
                                            json_data[model_id]["model_name"], hsp.bits, type(hsp.bits), type(true_pass_evalue)))

                                    blastResults = self.results(
                                        blastResults, blast_record.query + " | QUERY: " + alignment.hit_def, perfect, strict, loose, self.include_nudge)

            return blastResults
