"""Classes for DNA analysis (modified for Exercise 7)."""

from dna_functions import *

class Region(object):
    def __init__(self, dna, start, end):
        self._region = dna[start:end]

    def get_region(self):
        return self._region

    def __len__(self):
        return len(self._region)

    def __eq__(self, other):
        """Check if two Region instances are equal."""
        return self._region == other._region

    def __add__(self, other):
        """Add Region instances: self + other"""
        return self._region + other._region

    def __iadd__(self, other):
        """Increment Region instance: self += other"""
        self._region += other._region
        return self


class Gene(object):
    def __init__(self, dna=None, exon_regions=None):
        """
        Exercise 7:
        If called as Gene() with no arguments, generate random DNA sequence
        using generate_string from dna_functions.
        """

        # NEW: allow Gene() with no args
        if dna is None and exon_regions is None:
            # choose a reasonable default length (can be changed)
            dna = generate_string(60)
            exon_regions = None

        # --- original DNA handling ---
        if isinstance(dna, (list, tuple)) and \
           len(dna) == 2 and isinstance(dna[0], str) and \
           isinstance(dna[1], str):
            download(urlbase=dna[0], filename=dna[1])
            dna = read_dnafile(dna[1])
        elif isinstance(dna, str):
            pass
        else:
            raise TypeError(
                'dna=%s %s is not string or (urlbase,filename) tuple'
                % (dna, type(dna))
            )

        self._dna = dna

        er = exon_regions
        if er is None:
            self._exons = None
            self._introns = None
        else:
            if isinstance(er, (list, tuple)) and \
               len(er) == 2 and isinstance(er[0], str) and \
               isinstance(er[1], str):
                download(urlbase=er[0], filename=er[1])
                exon_regions = read_exon_regions(er[1])
            elif isinstance(er, (list, tuple)) and \
                 isinstance(er[0], (list, tuple)) and \
                 isinstance(er[0][0], int) and \
                 isinstance(er[0][1], int):
                pass
            else:
                raise TypeError(
                    'exon_regions=%s %s is not list of (int,int) '
                    'or (urlbase,filename) tuple' % (er, type(er))  # FIX: er not era
                )

            self._exon_regions = exon_regions

            # store exon Region objects
            self._exons = []
            for start, end in exon_regions:
                self._exons.append(Region(dna, start, end))

            # compute introns (regions between exons)
            self._introns = []
            prev_end = 0
            for start, end in exon_regions:
                if start - prev_end > 0:
                    self._introns.append(Region(dna, prev_end, start))
                prev_end = end
            if len(dna) - end > 0:
                self._introns.append(Region(dna, end, len(dna)))

    def get_dna(self):
        return self._dna

    def get_exons(self):
        return self._exons

    def get_introns(self):
        return self._introns