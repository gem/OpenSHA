// This class was generated by the JAXRPC SI, do not edit.
// Contents subject to change without notice.
// JAX-RPC Standard Implementation (1.1, build EA-R39)

package junk.webservices.client;

import com.sun.xml.rpc.encoding.DeserializationException;
import com.sun.xml.rpc.encoding.SOAPInstanceBuilder;
import com.sun.xml.rpc.util.exception.LocalizableExceptionAdapter;

public class GMT_WebServiceAPI_runGMT_Script_ResponseStruct_SOAPBuilder implements SOAPInstanceBuilder {
    private junk.webservices.client.GMT_WebServiceAPI_runGMT_Script_ResponseStruct _instance;
    private java.lang.String result;
    private static final int myRESULT_INDEX = 0;

    public GMT_WebServiceAPI_runGMT_Script_ResponseStruct_SOAPBuilder() {
    }

    public void setResult(java.lang.String result) {
        this.result = result;
    }

    public int memberGateType(int memberIndex) {
        switch (memberIndex) {
            case myRESULT_INDEX:
                return GATES_INITIALIZATION | REQUIRES_CREATION;
            default:
                throw new IllegalArgumentException();
        }
    }

    public void construct() {
    }

    public void setMember(int index, Object memberValue) {
        try {
            switch(index) {
                case myRESULT_INDEX:
                    _instance.setResult((java.lang.String)memberValue);
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
        catch (RuntimeException e) {
            throw e;
        }
        catch (Exception e) {
            throw new DeserializationException(new LocalizableExceptionAdapter(e));
        }
    }

    public void initialize() {
    }

    public void setInstance(Object instance) {
        _instance = (junk.webservices.client.GMT_WebServiceAPI_runGMT_Script_ResponseStruct)instance;
    }

    public Object getInstance() {
        return _instance;
    }
}