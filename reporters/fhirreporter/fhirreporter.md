# fhirreporter

FHIR® – Fast Healthcare Interoperability Resources (hl7.org/fhir) – is a next generation standards framework created by HL7. FHIR combines the best features of HL7's v2 icon, HL7 v3 icon and CDA icon product lines while leveraging the latest web standards and applying a tight focus on implementability.

FHIR solutions are built from a set of modular components called "Resources". These resources can easily be assembled into working systems that solve real-world clinical and administrative problems at a fraction of the price of existing alternatives. FHIR is suitable for use in a wide variety of contexts – mobile phone apps, cloud communications, EHR-based data sharing, server communication in large institutional healthcare providers, and much more.

This module produces a FHIR-format JSON files with variants' molecular consequences. For each sample in input, one JSON file (extension .fhir.json) will be created.

This module follows FHIR Genomics Operations.

Currently, molecluar consequence on Ensembl and RefSeq transcripts are returned. We are working with the CodeX FHIR Accelerator to incorporate more annotation data in FHIR.
