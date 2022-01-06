    Info<< "Create fluid mesh for time = "
        << runTime.timeName() << nl << endl;

const wordList fluidNNames
        (
	         rp.found("fluid") ? rp["fluid"] : wordList(0)
	);



    autoPtr<dynamicFvMesh> meshPtr		
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                "fluid",
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    dynamicFvMesh& mesh = meshPtr();		
