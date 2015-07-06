import re, sys, os


#TAXON_DEFINITION_BLOCK, SEQUENCE_DEFINITION_BLOCK, FILENAME

strict_constant = '''
<?xml version="1.0" standalone="yes"?>
<beast>
	<taxa id="taxa">
TAXON_DEFINITION_BLOCK
</taxa>

<alignment id="alignment" dataType="nucleotide">
SEQUENCE_DEFINITION_BLOCK
</alignment>

<patterns id="patterns" from="1" strip="false">
<alignment idref="alignment"/>
</patterns>

<constantSize id="constant" units="years">
<populationSize>
<parameter id="constant.popSize" value="150.0" lower="0.0"/>
</populationSize>
</constantSize>

<coalescentSimulator id="startingTree">
<taxa idref="taxa"/>
<constantSize idref="constant"/>
</coalescentSimulator>


	<treeModel id="treeModel">
		<coalescentTree idref="startingTree"/>
		<rootHeight>
			<parameter id="treeModel.rootHeight"/>
		</rootHeight>
		<nodeHeights internalNodes="true">
			<parameter id="treeModel.internalNodeHeights"/>
		</nodeHeights>
		<nodeHeights internalNodes="true" rootNode="true">
			<parameter id="treeModel.allInternalNodeHeights"/>
		</nodeHeights>
	</treeModel>

	<!-- Generate a coalescent likelihood                                        -->
	<coalescentLikelihood id="coalescent">
		<model>
			<constantSize idref="constant"/>
		</model>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</coalescentLikelihood>

	<!-- The strict clock (Uniform rates across branches)                        -->
	<strictClockBranchRates id="branchRates">
		<rate>
			<parameter id="clock.rate" value="1.0" lower="0.0"/>
		</rate>
	</strictClockBranchRates>

	<gtrModel id="gtr">
		<frequencies>
			<frequencyModel dataType="nucleotide">
				<frequencies>
					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>
				</frequencies>
			</frequencyModel>
		</frequencies>
		<rateAC>
			<parameter id="ac" value="1.0" lower="0.0"/>
		</rateAC>
		<rateAG>
			<parameter id="ag" value="1.0" lower="0.0"/>
		</rateAG>
		<rateAT>
			<parameter id="at" value="1.0" lower="0.0"/>
		</rateAT>
		<rateCG>
			<parameter id="cg" value="1.0" lower="0.0"/>
		</rateCG>
		<rateGT>
			<parameter id="gt" value="1.0" lower="0.0"/>
		</rateGT>
	</gtrModel>

	<siteModel id="siteModel">
		<substitutionModel>
			<gtrModel idref="gtr"/>
		</substitutionModel>
		<gammaShape gammaCategories="4">
			<parameter id="alpha" value="0.5" lower="0.0"/>
		</gammaShape>
	</siteModel>

	<treeLikelihood id="treeLikelihood" useAmbiguities="false">
		<patterns idref="patterns"/>
		<treeModel idref="treeModel"/>
		<siteModel idref="siteModel"/>
		<strictClockBranchRates idref="branchRates"/>
	</treeLikelihood>

	<operators id="operators" optimizationSchedule="default">
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="ac"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="ag"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="at"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="cg"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="gt"/>
		</scaleOperator>
		<deltaExchange delta="0.01" weight="0.1">
			<parameter idref="frequencies"/>
		</deltaExchange>
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="alpha"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="clock.rate"/>
		</scaleOperator>
		<subtreeSlide size="15.0" gaussian="true" weight="15">
			<treeModel idref="treeModel"/>
		</subtreeSlide>
		<narrowExchange weight="15">
			<treeModel idref="treeModel"/>
		</narrowExchange>
		<wideExchange weight="3">
			<treeModel idref="treeModel"/>
		</wideExchange>
		<wilsonBalding weight="3">
			<treeModel idref="treeModel"/>
		</wilsonBalding>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="treeModel.rootHeight"/>
		</scaleOperator>
		<uniformOperator weight="30">
			<parameter idref="treeModel.internalNodeHeights"/>
		</uniformOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="constant.popSize"/>
		</scaleOperator>
		<upDownOperator scaleFactor="0.75" weight="3">
			<up>
				<parameter idref="clock.rate"/>
			</up>
			<down>
				<parameter idref="treeModel.allInternalNodeHeights"/>
			</down>
		</upDownOperator>
	</operators>

	<mcmc id="mcmc" chainLength="100000000" autoOptimize="true">
		<posterior id="posterior">
			<prior id="prior">
				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
					<parameter idref="ac"/>
				</gammaPrior>
				<gammaPrior shape="0.05" scale="20.0" offset="0.0">
					<parameter idref="ag"/>
				</gammaPrior>
				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
					<parameter idref="at"/>
				</gammaPrior>
				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
					<parameter idref="cg"/>
				</gammaPrior>
				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
					<parameter idref="gt"/>
				</gammaPrior>
				<uniformPrior lower="0.0" upper="1.0">
					<parameter idref="frequencies"/>
				</uniformPrior>
				<exponentialPrior mean="0.5" offset="0.0">
					<parameter idref="alpha"/>
				</exponentialPrior>
				<uniformPrior lower="0.0" upper="1.0">
					<parameter idref="clock.rate"/>
				</uniformPrior>
				<oneOnXPrior>
					<parameter idref="constant.popSize"/>
				</oneOnXPrior>
<coalescentLikelihood idref="coalescent"/>
</prior>
<likelihood id="likelihood">
<treeLikelihood idref="treeLikelihood"/>
</likelihood>
</posterior>
<operators idref="operators"/>

<log id="screenLog" logEvery="5000">
<column label="Posterior" dp="4" width="12">
<posterior idref="posterior"/>
</column>
<column label="Prior" dp="4" width="12">
<prior idref="prior"/>
</column>
<column label="Likelihood" dp="4" width="12">
<likelihood idref="likelihood"/>
</column>
<column label="rootHeight" sf="6" width="12">
<parameter idref="treeModel.rootHeight"/>
</column>
<column label="clock.rate" sf="6" width="12">
<parameter idref="clock.rate"/>
</column>
</log>

<log id="fileLog" logEvery="5000" fileName="FILENAME_strict_constant.log" overwrite="false">
<posterior idref="posterior"/>
<prior idref="prior"/>
<likelihood idref="likelihood"/>
<parameter idref="treeModel.rootHeight"/>
<parameter idref="constant.popSize"/>
<parameter idref="ac"/>
<parameter idref="ag"/>
<parameter idref="at"/>
<parameter idref="cg"/>
<parameter idref="gt"/>
<parameter idref="frequencies"/>
<parameter idref="alpha"/>
<parameter idref="clock.rate"/>
<treeLikelihood idref="treeLikelihood"/>
<coalescentLikelihood idref="coalescent"/>
</log>

<logTree id="treeFileLog" logEvery="5000" nexusFormat="true" fileName="FILENAME_strict_constant.trees" sortTranslationTable="true">
<treeModel idref="treeModel"/>
<trait name="rate" tag="rate">
<strictClockBranchRates idref="branchRates"/>
</trait>
<posterior idref="posterior"/>
</logTree>
</mcmc>

<marginalLikelihoodEstimator chainLength="10000000" pathSteps="100" pathScheme="betaquantile" alpha="0.3">
<samplers>
<mcmc idref="mcmc"/>
</samplers>
<pathLikelihood id="pathLikelihood">
<source>
<posterior idref="posterior"/>
</source>
<destination>
<prior idref="prior"/>
</destination>
</pathLikelihood>
<log id="MLELog" logEvery="1000" fileName="FILENAME_strict_constant.mle.log">
<pathLikelihood idref="pathLikelihood"/>
</log>
</marginalLikelihoodEstimator>

<pathSamplingAnalysis fileName="FILENAME_strict_constant.mle.log">
<likelihoodColumn name="pathLikelihood.delta"/>
<thetaColumn name="pathLikelihood.theta"/>
</pathSamplingAnalysis>

<steppingStoneSamplingAnalysis fileName="FILENAME_strict_constant.mle.log">
<likelihoodColumn name="pathLikelihood.delta"/>
<thetaColumn name="pathLikelihood.theta"/>
</steppingStoneSamplingAnalysis>

<report>
<property name="timer">
<mcmc idref="mcmc"/>
</property>
</report>
</beast>


'''

#Note that TAXONNAME CAN BE SPLIT TO FIND THE DATE
taxon_block = '''
<taxon id="TAXONNAME">
<date value="DATE" direction="forwards" units="years"/>
</taxon>
'''

sequence_block = '''
<sequence>
<taxon idref="TAXONNAME"/>
SEQUENCE_DATA
</sequence>
'''
