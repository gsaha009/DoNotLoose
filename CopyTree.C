void CopyTree(TString fname)
{
   TString file = fname+"_FTree.root";
   std::cout<<"File: "<<file<<"\n";
   TFile *fin = TFile::Open(file);
   TTree *tin = dynamic_cast<TTree*>(fin->Get("RTree"));
   const auto nentries = tin->GetEntries();
   auto nTrain = static_cast<int>(nentries*0.70);
   std::cout<<"nEntries: "<<nentries<<"\t"<<"nTrain: "<<nTrain<<"\t"<<"nTest: "<<nentries-nTrain<<"\t";
   // Create a new file + a clone of old tree in new file
   TFile *outFile = new TFile(fname+".root", "recreate");
   TTree *train = (TTree*)tin->CloneTree(0);
   TTree *test  = (TTree*)tin->CloneTree(0);
   
   for (auto i : ROOT::TSeqI(nentries)) {
     tin->GetEntry(i);
     if (i <= nTrain) train->Fill();
     else test->Fill();
   }
   train->SetName("train");
   test ->SetName("test");
   outFile->Write();
}
