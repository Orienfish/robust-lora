/*
 * This program creates a simple network which uses an ADR algorithm to set up
 * the Spreading Factors of the devices in the Network.
 */

#include "ns3/point-to-point-module.h"
#include "ns3/forwarder-helper.h"
#include "ns3/network-server-helper.h"
#include "ns3/lora-channel.h"
#include "ns3/mobility-helper.h"
#include "ns3/lora-phy-helper.h"
#include "ns3/lorawan-mac-helper.h"
#include "ns3/lora-helper.h"
#include "ns3/gateway-lora-phy.h"
#include "ns3/periodic-sender.h"
#include "ns3/periodic-sender-helper.h"
#include "ns3/log.h"
#include "ns3/string.h"
#include "ns3/command-line.h"
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/lora-device-address-generator.h"
#include "ns3/random-variable-stream.h"
#include "ns3/config.h"
#include "ns3/rectangle.h"
#include "ns3/basic-energy-source-helper.h"
#include "ns3/lora-radio-energy-model-helper.h"
#include <numeric>

using namespace ns3;
using namespace lorawan;

NS_LOG_COMPONENT_DEFINE ("AdrExample");

// Trace sources that are called when a node changes its DR or TX power
void OnDataRateChange (uint8_t oldDr, uint8_t newDr)
{
  NS_LOG_DEBUG ("DR" << unsigned(oldDr) << " -> DR" << unsigned(newDr));
}

void OnTxPowerChange (double oldTxPower, double newTxPower)
{
  NS_LOG_DEBUG (oldTxPower << " dBm -> " << newTxPower << " dBm");
}

std::vector < std::string > split(std::string const & str, const char delim);
void EnablePeriodicDeviceStatusPrinting (NodeContainer endDevices,
                                         NodeContainer gateways,
                                         EnergySourceContainer sources,
                                         std::string filename,
                                         Time interval);
void DoPrintDeviceStatus (NodeContainer endDevices,
                          NodeContainer gateways,
                          EnergySourceContainer sources,
                          std::string filename);
double CalObjectiveValue (NodeContainer endDevices,
                          NodeContainer gateways, int M, 
                          EnergySourceContainer sources,
                          LoraPacketTracker& tracker, Time startTime, Time stopTime, double alpha);
double CalEnergyEfficiency (NodeContainer endDevices,
                            EnergySourceContainer sources,
                            std::string filename);
std::vector<double> CalEnergyEfficiency (NodeContainer endDevices,
                            EnergySourceContainer sources,
                            LoraPacketTracker& tracker, Time startTime, Time stopTime);


int main (int argc, char *argv[])
{

  bool verbose = false;
  bool adrEnabled = false;
  bool initializeSF = true;
  int nDevices = 1;
  int nPeriods = 1;
  double mobileNodeProbability = 0;
  double sideLength = 10000;
  int gatewayDistance = 5000;
  double maxRandomLoss = 10;
  int maxTransmissions = 8;
  std::string adrType = "ns3::AdrComponent";

  CommandLine cmd;
  cmd.AddValue ("verbose", "Whether to print output or not", verbose);
  cmd.AddValue ("MultipleGwCombiningMethod",
                "ns3::AdrComponent::MultipleGwCombiningMethod");
  cmd.AddValue ("MultiplePacketsCombiningMethod",
                "ns3::AdrComponent::MultiplePacketsCombiningMethod");
  cmd.AddValue ("HistoryRange", "ns3::AdrComponent::HistoryRange");
  cmd.AddValue ("MType", "ns3::EndDeviceLorawanMac::MType");
  cmd.AddValue ("EDDRAdaptation", "ns3::EndDeviceLorawanMac::EnableEDDataRateAdaptation");
  cmd.AddValue ("ChangeTransmissionPower",
                "ns3::AdrComponent::ChangeTransmissionPower");
  cmd.AddValue ("AdrEnabled", "Whether to enable ADR", adrEnabled);
  cmd.AddValue ("nDevices", "Number of devices to simulate", nDevices);
  cmd.AddValue ("PeriodsToSimulate", "Number of periods to simulate", nPeriods);
  cmd.AddValue ("sideLength",
                "Length of the side of the rectangle nodes will be placed in",
                sideLength);
  cmd.AddValue ("maxRandomLoss",
                "Maximum amount in dB of the random loss component",
                maxRandomLoss);
  cmd.AddValue ("gatewayDistance",
                "Distance between gateways",
                gatewayDistance);
  cmd.AddValue ("initializeSF",
                "Whether to initialize the SFs",
                initializeSF);
  cmd.AddValue ("MaxTransmissions",
                "Maximum number of retransmissions on end devices",
                maxTransmissions);
  cmd.Parse (argc, argv);


  /////////////
  // Logging //
  /////////////

  LogComponentEnable ("AdrExample", LOG_LEVEL_ALL);
  LogComponentEnable ("SimpleEndDeviceLoraPhy", LOG_LEVEL_ALL);
  LogComponentEnable ("LoraChannel", LOG_LEVEL_ALL); 
  LogComponentEnable ("PropagationLossModel", LOG_LEVEL_ALL);
  LogComponentEnable ("LoraPacketTracker", LOG_LEVEL_ALL);
  // LogComponentEnable ("NetworkServer", LOG_LEVEL_ALL);
  // LogComponentEnable ("NetworkController", LOG_LEVEL_ALL);
  // LogComponentEnable ("NetworkScheduler", LOG_LEVEL_ALL);
  // LogComponentEnable ("NetworkStatus", LOG_LEVEL_ALL);
  // LogComponentEnable ("EndDeviceStatus", LOG_LEVEL_ALL);
  // LogComponentEnable ("AdrComponent", LOG_LEVEL_ALL);
  // LogComponentEnable ("ClassAEndDeviceLorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable ("LogicalLoraChannelHelper", LOG_LEVEL_ALL);
  // LogComponentEnable ("MacCommand", LOG_LEVEL_ALL);
  // LogComponentEnable ("AdrExploraSf", LOG_LEVEL_ALL);
  // LogComponentEnable ("AdrExploraAt", LOG_LEVEL_ALL);
  // LogComponentEnable ("EndDeviceLorawanMac", LOG_LEVEL_ALL);
  LogComponentEnableAll (LOG_PREFIX_FUNC);
  LogComponentEnableAll (LOG_PREFIX_NODE);
  LogComponentEnableAll (LOG_PREFIX_TIME);

  // Set the EDs to require Data Rate control from the NS
  Config::SetDefault ("ns3::EndDeviceLorawanMac::DRControl", BooleanValue (true));



  //////////////////////////////////////
  // Create a simple wireless channel //
  //////////////////////////////////////

  Ptr<LogDistancePropagationLossModel> loss = CreateObject<LogDistancePropagationLossModel> ();
  loss->SetPathLossExponent (2.1);
  loss->SetReference (1000, 130);

  Ptr<UniformRandomVariable> x = CreateObject<UniformRandomVariable> ();
  x->SetAttribute ("Min", DoubleValue (0.0));
  x->SetAttribute ("Max", DoubleValue (maxRandomLoss));

  //Ptr<RandomPropagationLossModel> randomLoss = CreateObject<RandomPropagationLossModel> ();
  //randomLoss->SetAttribute ("Variable", PointerValue (x));

  //loss->SetNext (randomLoss);

  Ptr<PropagationDelayModel> delay = CreateObject<ConstantSpeedPropagationDelayModel> ();

  Ptr<LoraChannel> channel = CreateObject<LoraChannel> (loss, delay);



  ////////////////
  // Create EDs //
  ////////////////

  NodeContainer endDevices;
  endDevices.Create (nDevices);

  // End Device mobility
  MobilityHelper mobilityEd;
  Ptr<ListPositionAllocator> positionAllocEd = CreateObject<ListPositionAllocator> ();
  positionAllocEd->Add (Vector (0.0, 100000.0, 0.0));
  mobilityEd.SetPositionAllocator (positionAllocEd);
  mobilityEd.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  mobilityEd.Install (endDevices);

  // No exisiting end devices location file, then randomly generate locations and save it to the file
  /*std::string file_name = "EdLocation_" + std::to_string(nDevices) + ".txt";
  std::ifstream EdLocationFile(file_name);
  if (EdLocationFile.fail()) // no such file exists
  {
    NS_LOG_DEBUG ("Randomly generate ed device locations and save to file.");
    mobilityEd.SetPositionAllocator ("ns3::RandomRectanglePositionAllocator",
                                      "X", PointerValue (CreateObjectWithAttributes<UniformRandomVariable>
                                                         ("Min", DoubleValue(-sideLength),
                                                          "Max", DoubleValue(sideLength))),
                                      "Y", PointerValue (CreateObjectWithAttributes<UniformRandomVariable>
                                                         ("Min", DoubleValue(-sideLength),
                                                          "Max", DoubleValue(sideLength))));
    // Install mobility model on fixed nodes and output location
    mobilityEd.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
    std::ofstream OutputFile(file_name);
    for (int i = 0; i < nDevices; ++i) {
      mobilityEd.Install (endDevices.Get (i));
      Ptr<MobilityModel> position = endDevices.Get (i) -> GetObject<MobilityModel> ();
      Vector pos = position->GetPosition ();
      OutputFile << pos.x << " " << pos.y << std::endl;
    }
    OutputFile.close();   
  }
  else // read from file
  {
    NS_LOG_DEBUG ("Read from existing ed device location file.");
    std::string line;
    Ptr<ListPositionAllocator> positionAllocEd = CreateObject<ListPositionAllocator> ();
    while (std::getline(EdLocationFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > coordinates = split(line, ' ');
            double x = atof(coordinates.at(0).c_str());
            double y = atof(coordinates.at(1).c_str());
            positionAllocEd->Add (Vector (x, y, 0.0) );
        }
    }
    mobilityEd.SetPositionAllocator (positionAllocEd);
    mobilityEd.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
    mobilityEd.Install (endDevices);
  }*/

  // Install mobility model on mobile nodes
  //mobilityEd.SetMobilityModel ("ns3::RandomWalk2dMobilityModel",
  //                             "Bounds", RectangleValue (Rectangle (-sideLength, sideLength,
  //                                                                  -sideLength, sideLength)),
  //                             "Distance", DoubleValue (1000),
  //                             "Speed", PointerValue (CreateObjectWithAttributes<UniformRandomVariable>
  //                                                    ("Min", DoubleValue(minSpeed),
  //                                                     "Max", DoubleValue(maxSpeed))));
  //for (int i = fixedPositionNodes; i < (int) endDevices.GetN (); ++i)
  //  {
  //    mobilityEd.Install (endDevices.Get (i));
  //  }



  ////////////////
  // Create GWs //
  ////////////////

  NodeContainer gateways;
  int nGateways = 1;
  gateways.Create (nGateways);

  MobilityHelper mobilityGw;
  Ptr<ListPositionAllocator> positionAllocGw = CreateObject<ListPositionAllocator> ();
  positionAllocGw->Add (Vector (0.0, 0.0, 15.0));
  //positionAllocGw->Add (Vector (-5000.0, -5000.0, 15.0));
  //positionAllocGw->Add (Vector (-5000.0, 5000.0, 15.0));
  //positionAllocGw->Add (Vector (5000.0, -5000.0, 15.0));
  //positionAllocGw->Add (Vector (5000.0, 5000.0, 15.0));
  mobilityGw.SetPositionAllocator (positionAllocGw);
  mobilityGw.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  //Ptr<HexGridPositionAllocator> hexAllocator = CreateObject<HexGridPositionAllocator> (gatewayDistance / 2);
  //mobilityGw.SetPositionAllocator (hexAllocator);
  //mobilityGw.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  mobilityGw.Install (gateways);



  /////////////////////////////
  // Create the LoRa Helpers //
  /////////////////////////////

  LoraPhyHelper phyHelper = LoraPhyHelper ();
  phyHelper.SetChannel (channel);

  // Create the LorawanMacHelper
  LorawanMacHelper macHelper = LorawanMacHelper ();

  // Create the LoraHelper
  LoraHelper helper = LoraHelper ();
  helper.EnablePacketTracking ();

  // Create the LoraNetDevices of the gateways
  phyHelper.SetDeviceType (LoraPhyHelper::GW);
  macHelper.SetDeviceType (LorawanMacHelper::GW);
  helper.Install (phyHelper, macHelper, gateways);

  // Create a LoraDeviceAddressGenerator
  uint8_t nwkId = 54;
  uint32_t nwkAddr = 1864;
  Ptr<LoraDeviceAddressGenerator> addrGen = CreateObject<LoraDeviceAddressGenerator> (nwkId,nwkAddr);

  // Create the LoraNetDevices of the end devices
  phyHelper.SetDeviceType (LoraPhyHelper::ED);
  macHelper.SetDeviceType (LorawanMacHelper::ED_A);
  macHelper.SetAddressGenerator (addrGen);
  macHelper.SetRegion (LorawanMacHelper::EU);
  NetDeviceContainer endDevicesNetDevices = helper.Install (phyHelper, macHelper, endDevices);



  /////////////////////////////////
  // Install applications in EDs //
  /////////////////////////////////

  int appPeriodSeconds = 1200;      // One packet every 20 minutes
  PeriodicSenderHelper appHelper = PeriodicSenderHelper ();
  appHelper.SetPeriod (Seconds (appPeriodSeconds));
  ApplicationContainer appContainer = appHelper.Install (endDevices);

  // Do not set spreading factors up: we will wait for the NS to do this
  if (initializeSF) 
  {
    macHelper.SetSpreadingFactorsUp (endDevices, gateways, channel);
  }



  ///////////////
  // Create NS //
  ///////////////

  NodeContainer networkServers;
  networkServers.Create (1);

  // Install the NetworkServer application on the network server
  NetworkServerHelper networkServerHelper;
  networkServerHelper.SetGateways (gateways);
  networkServerHelper.SetEndDevices (endDevices);
  networkServerHelper.EnableAdr (adrEnabled);
  networkServerHelper.SetAdr (adrType);
  networkServerHelper.Install (networkServers);

  // Install the Forwarder application on the gateways
  ForwarderHelper forwarderHelper;
  forwarderHelper.Install (gateways);

  //////////////////////////
  // Install Energy Model //
  //////////////////////////

  BasicEnergySourceHelper basicSourceHelper;
  LoraRadioEnergyModelHelper radioEnergyHelper;

  // configure energy source
  basicSourceHelper.Set ("BasicEnergySourceInitialEnergyJ", DoubleValue (10000)); // Energy in J
  basicSourceHelper.Set ("BasicEnergySupplyVoltageV", DoubleValue (3.3));

  radioEnergyHelper.Set ("StandbyCurrentA", DoubleValue (0.0014));
  radioEnergyHelper.Set ("TxCurrentA", DoubleValue (0.028));
  radioEnergyHelper.Set ("SleepCurrentA", DoubleValue (0.0000015));
  radioEnergyHelper.Set ("RxCurrentA", DoubleValue (0.0112));

  radioEnergyHelper.SetTxCurrentModel ("ns3::ConstantLoraTxCurrentModel",
                                       "TxCurrent", DoubleValue (0.028));

  // install source on EDs' nodes
  EnergySourceContainer sources = basicSourceHelper.Install (endDevices);
  Names::Add ("/Names/EnergySource", sources.Get (0));

  // install device model
  DeviceEnergyModelContainer deviceModels = radioEnergyHelper.Install
      (endDevicesNetDevices, sources);



  /////////////
  // Tracing //
  /////////////
  // Connect our traces
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/0/$ns3::LoraNetDevice/Mac/$ns3::EndDeviceLorawanMac/TxPower",
                                 MakeCallback (&OnTxPowerChange));
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/0/$ns3::LoraNetDevice/Mac/$ns3::EndDeviceLorawanMac/DataRate",
                                 MakeCallback (&OnDataRateChange));

  // Activate printing of ED MAC parameters
  Time stateSamplePeriod = Seconds (1200);
  EnablePeriodicDeviceStatusPrinting (endDevices, gateways, sources, "nodeData.txt", stateSamplePeriod);
  helper.EnablePeriodicPhyPerformancePrinting (gateways, "phyPerformance.txt", stateSamplePeriod);
  helper.EnablePeriodicGlobalPerformancePrinting ("globalPerformance.txt", stateSamplePeriod);

  // Activate packet tracker
  LoraPacketTracker& tracker = helper.GetPacketTracker ();

  // Start simulation
  Time simulationTime = Seconds (1200 * nPeriods);
  Simulator::Stop (simulationTime);
  Simulator::Run ();
  Simulator::Destroy ();

  int M = 1;
  double alpha = 0.5;
  std::cout << tracker.CountMacPacketsGlobally (Seconds (0), simulationTime) << std::endl;
  std::cout << tracker.CountMacPacketsGloballyCpsr (Seconds (0), simulationTime) << std::endl;
  std::cout << tracker.PrintPhyPacketsPerGw (Seconds (0), simulationTime, nDevices) << std::endl;
  std::cout << CalObjectiveValue (endDevices, gateways, M, sources, tracker, Seconds (0), simulationTime, alpha) << std::endl;

  return 0;
}


// split implementation for reading from external text files
std::vector < std::string > split(std::string const & str, const char delim) 
{
    std::vector < std::string > result;

    std::stringstream ss(str);
    std::string s;

    while (std::getline(ss, s, delim)) {
        result.push_back(s);
    }

    return result;
}

// Schedule periodic end devices' status printing - including energy
void EnablePeriodicDeviceStatusPrinting (NodeContainer endDevices,
                                         NodeContainer gateways,
                                         EnergySourceContainer sources,
                                         std::string filename,
                                         Time interval)
{
  DoPrintDeviceStatus (endDevices, gateways, sources, filename);

  // Schedule periodic printing
  Simulator::Schedule (interval,
                       &EnablePeriodicDeviceStatusPrinting,
                       endDevices, gateways, sources, filename, interval);
}

// Event function of printing end devices' status - including energy
void DoPrintDeviceStatus (NodeContainer endDevices,
                          NodeContainer gateways,
                          EnergySourceContainer sources,
                          std::string filename)
{
  const char * c = filename.c_str ();
  std::ofstream outputFile;
  if (Simulator::Now () == Seconds (0))
  {
    // Delete contents of the file as it is opened
    outputFile.open (c, std::ofstream::out | std::ofstream::trunc);
  }
  else
  {
    // Only append to the file
    outputFile.open (c, std::ofstream::out | std::ofstream::app);
  }

  Time currentTime = Simulator::Now();
  
  // print end devices' status
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
  {
    Ptr<Node> object = *j;
    int nodeId = object->GetId();

    // Get position
    Ptr<MobilityModel> position = object->GetObject<MobilityModel> ();
    NS_ASSERT (position != 0);
    Vector pos = position->GetPosition ();

    // Get netDevice
    Ptr<NetDevice> netDevice = object->GetDevice (0);
    Ptr<LoraNetDevice> loraNetDevice = netDevice->GetObject<LoraNetDevice> ();
    NS_ASSERT (loraNetDevice != 0);

    // Get mac parameters
    Ptr<ClassAEndDeviceLorawanMac> mac = loraNetDevice->GetMac ()->GetObject<ClassAEndDeviceLorawanMac> ();
    int dr = int(mac->GetDataRate ());
    double txPower = mac->GetTransmissionPower ();

    // Get energy
    Ptr<EnergySource> SourcePtr = sources.Get(nodeId);
    NS_ASSERT (SourcePtr != 0);
    Ptr<DeviceEnergyModel> LoRaRadioModelPtr = SourcePtr->FindDeviceEnergyModels("ns3::LoraRadioEnergyModel").Get(0);
    NS_ASSERT (LoRaRadioModelPtr != 0);
    double energy_consumption = LoRaRadioModelPtr->GetTotalEnergyConsumption();
    
    outputFile << currentTime.GetSeconds () << " "
               << object->GetId () <<  " "
               << pos.x << " " << pos.y << " " << pos.z << " "
               << dr << " " << unsigned(txPower) << " "
               << energy_consumption << std::endl;
  }

  // print gateways' status
  for (NodeContainer::Iterator j = gateways.Begin (); j != gateways.End (); ++j)
  {
    Ptr<Node> object = *j;
    Ptr<MobilityModel> position = object->GetObject<MobilityModel> ();
    Vector pos = position->GetPosition ();
    outputFile << currentTime.GetSeconds () << " "
               << object->GetId () <<  " "
               << pos.x << " " << pos.y << " " << pos.z << " "
               << "-1 -1 -1" << std::endl;
  }
  outputFile.close ();
}

// Calculate the objective value in the ICIOT paper
double CalObjectiveValue (NodeContainer endDevices,
                          NodeContainer gateways, int M, 
                          EnergySourceContainer sources,
                          LoraPacketTracker& tracker, Time startTime, Time stopTime, double alpha)
{
  std::vector<double> vecEE = CalEnergyEfficiency(endDevices, sources, tracker, startTime, stopTime);
  int edCnt = endDevices.GetN();
  double EE = std::accumulate(vecEE.begin(), vecEE.end(), 0) / edCnt;

  int gwCnt = gateways.GetN();
  return EE - alpha / M * gwCnt;
}

// Calculate the energy efficiency across all end devices
std::vector<double> CalEnergyEfficiency (NodeContainer endDevices,
                            EnergySourceContainer sources,
                            LoraPacketTracker& tracker, Time startTime, Time stopTime)
{
  std::vector<double> vecEE;

  // iterative through all end devices
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
  {
    Ptr<Node> object = *j;
    int edId = object->GetId();
    NS_LOG_DEBUG ("Dealing with node ID " << edId);

    // Get packet delivery ratio
    std::vector<int> packetEd = tracker.CountPhyPacketsPerEd (startTime, stopTime, edId);
    double pdr = packetEd.at (2) / packetEd.at (1);
    NS_LOG_DEBUG ("Packet statistics: sent " << packetEd.at(1) << " received " << packetEd.at(2));
    NS_LOG_DEBUG ("Packet delivery ratio " << pdr);

    // Get energy
    Ptr<EnergySource> SourcePtr = sources.Get(edId);
    NS_ASSERT (SourcePtr != 0);
    Ptr<DeviceEnergyModel> LoRaRadioModelPtr = SourcePtr->FindDeviceEnergyModels("ns3::LoraRadioEnergyModel").Get(0);
    NS_ASSERT (LoRaRadioModelPtr != 0);
    double energyPerPkt = LoRaRadioModelPtr->GetTotalEnergyConsumption() / packetEd.at (1);
    NS_LOG_DEBUG ("Energy consumption per sent packet: " << energyPerPkt);

    // Push back energy efficiency
    double nodeEE = pdr / energyPerPkt;
    NS_LOG_DEBUG ("Energy efficiency: " << nodeEE);
    vecEE.push_back (nodeEE);
  }

  return vecEE;
}
