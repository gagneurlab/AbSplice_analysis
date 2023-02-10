from typing import List
from kipoiseq import Interval
from pyranges import PyRanges
import pandas as pd
from typing import List, Union, Iterable, Iterator


class Junction(Interval):

    @property
    def acceptor(self):
        return self.start if self.strand == '-' else self.end

    @property
    def donor(self):
        return self.end if self.strand == '-' else self.start

    def dinucleotide_region(self):
        return Interval(self.chrom, self.start, self.start + 2), \
            Interval(self.chrom, self.end - 2, self.end)

    def acceptor_region(self, overhang=(250, 250)):
        return Interval(self.chrom, self.acceptor,
                        self.acceptor, strand=self.strand) \
            .slop(upstream=overhang[0], downstream=overhang[1])

    def donor_region(self, overhang=(250, 250)):
        return Interval(self.chrom, self.donor,
                        self.donor, strand=self.strand) \
            .slop(upstream=overhang[0], downstream=overhang[1])


def get_splice_site_intervals(junction, overhang=(250, 250)):
    junction = Junction.from_str(junction) if type(
        junction) == str else junction

    acceptor = junction.acceptor_region(overhang=overhang)
    donor = junction.donor_region(overhang=overhang)
    return [acceptor, donor]

def get_unique_splice_site_intervals_in_event(event, overhang=(250, 250)):
    sites = list()
    for junction in event:
        sites.append(get_splice_site_intervals(junction, overhang))
    sites = [item for sublist in sites for item in sublist]
    sites = list(set(sites))
    return sites
    
def intervals_to_pyranges(intervals: List[Interval]) -> PyRanges:
    """
    Create pyrange object given list of intervals objects.
    Args:
      intervals: list of interval objects have CHROM, START, END, properties.
    """
    import pyranges
    df = pd.DataFrame([
        (
            i.chrom,
            i.start,
            i.end,
            i
        )
        for i in intervals
    ], columns=['Chromosome', 'Start', 'End', 'interval'])
    return pyranges.PyRanges(df)
    
    
    
class Event:

    def __init__(self, junctions: List):
        self.junctions = junctions

    @classmethod
    def from_str(cls, event_str: str):
        return cls([
            Junction.from_str(junc)
            for junc in event_str.split(';')
        ])


class EventPSI5(Event):

    def __init__(self, junctions: List):
        super().__init__(junctions)
        assert all(
            i.donor == self.junctions[0].donor
            for i in self.junctions
        )

    @property
    def donor_str(self):
        j = self.junctions[0]
        return f'{j.chrom}:{j.donor}:{j.strand}'


class EventPSI3(Event):

    def __init__(self, junctions: List):
        super().__init__(junctions)
        assert all(
            i.acceptor == self.junctions[0].acceptor
            for i in self.junctions
        )

    @property
    def acceptor_str(self):
        j = self.junctions[0]
        return f'{j.chrom}:{j.acceptor}:{j.strand}'
