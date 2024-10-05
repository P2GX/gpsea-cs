import typing

import hpotk

from gpsea.model import Patient
from gpsea.analysis.pscore import PhenotypeScorer


class Ankrd11PhenotypeScorer(PhenotypeScorer):
    """
    `Ankrd11PhenotypeScorer` computes the phenotype score
    as described in `Martinez-Cayuelas et al. <https://pubmed.ncbi.nlm.nih.gov/36446582>`_.
    """

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        self._hpo = hpo

        # severe and profound GDD
        self._gdd_tids = {
            'HP:0011344': 2, 'HP:0012736': 2,
            'HP:0011342': 1, 'HP:0011343': 1, 'HP:0001263': 1,
        }

       

    def _developmental_delay_score(
        self,
        observed_term_ids: typing.Iterable[str],
    ) -> float:
        """
        1. History of developmental delay:
        1.1. Language delay
        1.2. Motor delay 
        
        Args:
            observed_term_ids: terms observed in patient

        Returns: 1 point of GDD, motor delay, or language delay is present
        """
        # Check GDD terms with higher priority than ID terms.
        # Global developmental delay
        gdd = hpotk.TermId.from_curie("HP:0001263") # Global developmental delay
        md = hpotk.TermId.from_curie("HP:0001270") # Motor delay
        sld = hpotk.TermId.from_curie("HP:0000750") # Delayed speech and language development
        dev_delay_terms = [gdd, md, sld]
        for t in observed_term_ids:
            for ddt in dev_delay_terms:
                if t == gdd or self._hpo.graph.is_descendant_of(t, ddt):
                    return 1
        return 0
    
    def _id_asd_adhd_score(self,
        observed_term_ids: typing.Iterable[str],) -> float:
        """
        ID and/or ASD and/or ADHD
        """
         # mild, moderate, and unspecified 1D (borderline has 0.5)
        idd_tids = {
            'HP:0010864': 2, 'HP:0002187': 2,
            'HP:0001256': 1, 'HP:0002342': 1, 'HP:0001249': 1,
            'HP:0006889': 0.5,
        }
        adhd = hpotk.TermId.from_curie("HP:0007018") # Attention deficit hyperactivity disorder
        aut_behavior = hpotk.TermId.from_curie("HP:0000729") # Autistic behavior
        for t in observed_term_ids:
            if t in idd_tids:
                return 1
            if t == adhd:
                return 1
            if t == aut_behavior or self._hpo.graph.is_descendant_of(t, aut_behavior):
                return 1
        return 0
            
     

    def _term_or_descendant(
        self,
        target_tid: str,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        Args:
            target_tid: term of interest
            observed_term_ids: all terms observed in patient

        Returns:
            1 if the term or any descendant is present in the patient, otherwise 0
        """
        for term_id in observed_term_ids:
            if term_id == target_tid \
               or any(ancestor == target_tid for ancestor in self._hpo.graph.get_ancestors(term_id)):
                return 1
        
        return 0

    def _term_or_descendant_count(
        self,
        target_tid: str,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        Args:
            target_tid: term of interest
            observed_term_ids: all terms observed in patient

        Returns:
            the total count of the terms equal to or descending from the target_tid
        """
        total_count = 0
        for term_id in observed_term_ids:
            for desc_tid in self._hpo.graph.get_ancestors(term_id, include_source=True):
                if desc_tid.value == target_tid:
                    total_count += 1
        if total_count > 0:
            return 1
        else:
            return 0

    def _comorbidities(
        self,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        One point each for
        3.1. Hearing loss
        3.2. Recurrent otitis media
        3.3. Visual problems
        3.4. Congenital heart defects
        3.5. Seizures
        3.6. Feeding difficulties
        3.7. Cryptorchidism 
        
        Args:
            observed_term_ids: terms observed in patient

        Returns: an `int` (between 0 and 7)
        """
        hearing = hpotk.TermId.from_curie("HP:0000365") # Hearing impairment
        otitis  = hpotk.TermId.from_curie("HP:0000403") #Recurrent otitis media
        visual  = hpotk.TermId.from_curie("HP:0000505")  # Visual impairment 
        heart   = hpotk.TermId.from_curie("HP:0001627")  # Abnormal heart morphology, approximation
        seizure = hpotk.TermId.from_curie("HP:0001250") # Seizure  
        feeding = hpotk.TermId.from_curie("HP:0011968") # Feeding difficulties
        crypt = hpotk.TermId.from_curie("HP:0000028") # Cryptorchidism 
        total_count = 0
        for t in observed_term_ids:
            for tid in (hearing, otitis, visual, heart, seizure, feeding, crypt):
                total_count += self._term_or_descendant_count(tid, 1)
        return total_count



    def _phenotypic_score(
        self,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        4.1 Macrodontia and/or other dental anomalies
        4.2. Triangular face
        4.3. Long philtrum
        4.4. Characteristic nose (anteverted nares, and/or bulbous tip, and/or
        prominent nose)
        4.5. Characteristic eyebrows (wide and/or thick, and/or bushy, and/or
        synophrys)
        4.6. Characteristic ears (large and/or prominent and/or low-set)
        4.7. Hand anomalies (brachydactyly and/or clinodactyly)
        4.8. Postnatal short stature 
        One point is assigned for either the corresponding HPO terms or any of their descendents.
        
        Args:
            observed_term_ids:  terms observed in patient

        Returns:   Non-facial dysmorphism and congenital abnormalities score (between 0 and 2)

        """
        tooth = hpotk.TermId.from_curie("HP:0006482") # Abnormal dental morphology
        triangular = hpotk.TermId.from_curie("HP:0000325") #Triangular face 
        philtrum = hpotk.TermId.from_curie("HP:0000343") # Long philtrum 

        nose_d = {
            "HP:0000463":"Anteverted nares",
            "HP:0000414": "Bulbous nose",
            "HP:0000448":"Prominent nose"
        }

        ears_d = {
            "HP:0000400": "Macrotia",
            "HP:0000411": "Protruding ear",
            "HP:0000369": "Low-set ears"
        }

        hand = hpotk.TermId.from_curie("HP:0005922") # Abnormal hand morphology 
        short = hpotk.TermId.from_curie("HP:0004322") #  Short stature 

        items = [tooth, triangular, philtrum, hand, short]
        pheno_score = 0
        for t in observed_term_ids:
            for item in items:
                pheno_score += self._term_or_descendant_count(item)
            if t.value in ears_d:
                pheno_score += 1
            elif t.value in nose_d:
                pheno_score += 1
            
        return pheno_score

   

    def score(self, patient: Patient) -> float:
        """
        Calculate score based on list of strings with term identifiers or observed HPO terms.
        
        Args:
            patient: list of strings with term identifiers or observed HPO terms

        Returns: de Vries score between 0 and 10

        """
        observed_term_ids = tuple(tid.identifier.value for tid in patient.present_phenotypes())

        phenoscore = 0
        phenoscore += self._developmental_delay_score(observed_term_ids=observed_term_ids)
        phenoscore += self._id_asd_adhd_score(observed_term_ids=observed_term_ids)
        phenoscore += self._comorbidities(observed_term_ids=observed_term_ids)
        phenoscore += self._phenotypic_score(observed_term_ids=observed_term_ids)
        
        return phenoscore