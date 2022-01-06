    Info<< "Create solid mesh for time = "
        << runTime.timeName() << nl << endl;

const wordList solidNames
        (
	         rp.found("solid") ? rp["solid"] : wordList(0)
	);



    autoPtr<dynamicFvMesh> mesh2Ptr		//insert 2 for solid
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                "solid",
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    dynamicFvMesh& mesh2 = mesh2Ptr();		//insert 2 for solid
