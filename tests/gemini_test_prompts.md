# Gemini CLI Test Prompts for Rosetta MCP Server

## Setup Verification

The Rosetta MCP server has been successfully added to Gemini CLI configuration at `~/.gemini/settings.json`.

Configuration details:
- **Server name**: `rosetta`
- **Command**: `/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/env/bin/python`
- **Script**: `/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/src/server.py`

## Test Prompts for Gemini CLI

Use these prompts to test the Rosetta MCP integration:

### Basic Tool Discovery
```
What MCP tools are available from the rosetta server? List them with descriptions.
```

### Structure Validation
```
Use the rosetta server to validate this PDB structure: /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/examples/data/test_input.pdb
```

### Example Structures
```
What example structures are available for testing with the rosetta server?
```

### Job Management
```
List all current jobs using the rosetta server.
```

### Error Handling Test
```
Try to validate a non-existent file using rosetta: /fake/nonexistent.pdb
```

### Complete Workflow
```
Using the rosetta server:
1. Validate the structure /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/examples/data/test_input.pdb
2. If valid, refine it with 3 trajectories
3. Show me the results
```

### Advanced Testing
```
Submit a protein docking job using the rosetta server for the file /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/examples/data/test_docking.pdb and monitor its progress.
```

## Usage Instructions

To test with Gemini CLI:

1. **Ensure API access**: Make sure you have valid Gemini API credentials configured
2. **Start Gemini CLI**: Run `gemini` in your terminal
3. **Use prompts**: Copy and paste the test prompts above
4. **Interactive mode**: You can also use `gemini -i` for interactive mode

## Expected Behavior

- **Tool Discovery**: Should list all 14 Rosetta tools
- **Structure Validation**: Should return validation results with atom counts, chains, etc.
- **Error Handling**: Should return structured error messages for invalid inputs
- **Job Submission**: Should return job IDs for async operations
- **Job Monitoring**: Should allow checking status and retrieving results

## Troubleshooting

### Server Connection Issues
If Gemini CLI can't connect to the rosetta server:
1. Verify the Python environment path is correct
2. Test the server manually: `python /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/src/server.py`
3. Check Gemini settings: `cat ~/.gemini/settings.json`

### API Rate Limits
If you encounter "model overloaded" errors:
- Wait a few minutes and try again
- Use shorter, simpler prompts
- Consider using --model flag with different model variants

## Integration Status

✅ **Server Added**: Successfully added to Gemini CLI configuration
✅ **Connection Test**: MCP client initialization successful
⚠️ **API Testing**: Limited by current API availability
📝 **Manual Testing**: Required when API is available

## Notes

- The Gemini CLI will automatically discover and use Rosetta MCP tools
- Tool parameters are validated by the MCP protocol
- Results are returned in structured JSON format
- Long-running jobs can be monitored through job management tools